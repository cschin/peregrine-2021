// Peregrine Assembler and SHIMMER Genome Assembly Toolkit
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

const VERSION_STRING: &'static str = env!("VERSION_STRING");

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
//static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

use clap::clap_app;
use glob::glob;
use std::fs::File;
use std::fs::{create_dir_all, remove_file};
use std::io::{BufRead, Write};
use std::io::{BufReader, BufWriter, Result};
use std::path::Path;
use sysinfo::SystemExt;

mod utils;
use simple_logger::SimpleLogger;
use std::time::SystemTime;
use utils::build_idx;
use utils::build_sdb;
use utils::dp_graph;
use utils::graph;
use utils::layout;
use utils::ovlp;
use utils::ovlp_ec;
use utils::resolve::resolve_ht;
use utils::seqmap::dedup_target_seqs;
use utils::Parameters;
use utils::{getrusage, log_resource, rusage, MaybeUninit, RUSAGE_SELF};

fn cat_path(wd: &String, filename: &String) -> String {
    Path::new(&wd).join(&filename).to_string_lossy().to_string()
}

fn get_ovlps(
    input_reads: &String,
    work_dir: &String,
    prefix: &String,
    delete_input: bool,
    parameters: &Parameters,
    rdata: &mut rusage,
) -> Result<()> {
    log_resource(
        &format!("BGN: get_ovlp, input_reads: {}", input_reads),
        rdata,
    );
    let output_prefix = format!("{}/{}", &work_dir, &prefix);

    log_resource("BGN: bulding sequence database", rdata);
    let nbase = build_sdb::build(&input_reads, &output_prefix)?;
    log_resource("END: bulding sequence database", rdata);
    if delete_input {
        let seq_list_buf = BufReader::new(File::open(&input_reads).unwrap());
        for fastx_file in seq_list_buf.lines() {
            let fastx_file = fastx_file?;
            remove_file(&fastx_file)?;
        }
    }

    let seqdb = cat_path(&work_dir, &format!("{}.seqdb", &prefix));
    let seqidx = cat_path(&work_dir, &format!("{}.idx", &prefix));
    let shmmer_idx = cat_path(&work_dir, &format!("{}-shmr", &prefix));

    // step 2: build shimmer index
    log_resource("BGN: bulding shimmer", rdata);
    build_idx::build(&seqdb, &seqidx, &shmmer_idx, &parameters);
    log_resource("END: bulding shimmer", rdata);

    //step 3: overlap reads
    let ovlp_out = format!("{}/{}-ovlp", &work_dir, &prefix);
    log_resource("BGN: overlapping", rdata);
    let system = sysinfo::System::new_all();
    let free_mem = system.total_memory() - system.used_memory();
    if (free_mem as f64) < ((nbase >> 10) as f64 * 1.5) {
        log::warn!("free memory is less than 1.5 x (total number of bases)");
        log::warn!(
            "free memory = {}kb, total bases / 1024 = {}",
            free_mem,
            nbase >> 10
        );
    }
    ovlp::ovlp(&seqdb, &seqidx, &shmmer_idx, &ovlp_out, &parameters)?;
    log_resource("END: overlapping", rdata);

    log_resource(
        &format!("END: get_ovlp - input_reads: {}", input_reads),
        rdata,
    );
    Ok(())
}

fn main() -> Result<()> {
    let mut rdata: MaybeUninit<libc::rusage> = unsafe { MaybeUninit::uninit().assume_init() };
    let _res = unsafe { getrusage(RUSAGE_SELF, &mut rdata.assume_init_read()) };

    let matches = clap_app!(pg_asm =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler
pg_asm: the main workflow entry for end-to-end assembly from the reads
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")
        (@arg input_reads: +required "Path to a file that contains the list of reads in .fa .fa.gz .fastq or fastq.gz formats")
        (@arg work_dir: +required "The path to a work directory for intermediate files and the results")
        (@arg NTHREADS: +takes_value "Number of threads")
        (@arg NCHUNKS: +takes_value "Number of partition")
        (@arg w: -w +takes_value "Window size [default: 80]")
        (@arg k: -k +takes_value "Kmer size [default: 56]")
        (@arg r: -r +takes_value "Reduction factor [default: 6]")
        (@arg tol: -t --tol +takes_value "Alignment tolerance [default: 0.01]")
        (@arg layout_method: -l +takes_value "layout version [default: 2]")
        (@arg bestn: --bestn -b +takes_value "number of best overlaps for initial graph [default: 6]")
        (@arg keep: --keep "keep intermediate files")
        (@arg fast: --fast "run the assembler in the fast mode")
        (@arg no_resolve: --no_resolve "disable resolving repeats / dups at the end")
        (@arg min_ec_cov: -c --min_ec_cov +takes_value "Minimum error coverage [default: 1]")
        (@arg log: --log +takes_value "log level: DBBUG or INFO (default)")
    ).get_matches();

    let log_level = match matches.value_of("log").unwrap_or("INFO") {
        "DEBUG" => log::LevelFilter::Debug,
        _ => log::LevelFilter::Info,
    };

    SimpleLogger::new()
        .with_level(log_level)
        .with_utc_timestamps()
        .init()
        .unwrap();

    let input_reads = matches.value_of("input_reads").unwrap().to_string();
    let work_dir = matches.value_of("work_dir").unwrap().to_string();

    let prefix = "reads".to_string();

    let keep = matches.is_present("keep");
    let fastmode = matches.is_present("fast");
    let no_resolve = matches.is_present("no_resolve");

    let physical_cpus = num_cpus::get_physical();
    let nthreads = matches
        .value_of("NTHREADS")
        .unwrap_or(&physical_cpus.to_string())
        .parse::<u32>()
        .unwrap();

    let nchunks: u32;
    if matches.is_present("NCHUNKS") {
        nchunks = matches
            .value_of("NCHUNKS")
            .unwrap()
            .to_string()
            .parse::<u32>()
            .unwrap();
    } else {
        let physical_cpus = physical_cpus as u32;
        match physical_cpus {
            1..=5 => nchunks = 16 + physical_cpus * 4,
            6..=12 => nchunks = 2 + physical_cpus * 3,
            13..=19 => nchunks = physical_cpus * 2,
            _ => nchunks = physical_cpus,
        }
    }

    let wsize = matches
        .value_of("w")
        .unwrap_or("80")
        .parse::<u32>()
        .unwrap();

    let ksize = matches
        .value_of("k")
        .unwrap_or("56")
        .parse::<u32>()
        .unwrap();

    let rfactor = matches.value_of("r").unwrap_or("6").parse::<u32>().unwrap();

    let tol = matches
        .value_of("tol")
        .unwrap_or("0.01")
        .parse::<f64>()
        .unwrap();

    let min_ec_cov = matches
        .value_of("min_ec_cov")
        .unwrap_or("1")
        .parse::<u16>()
        .unwrap();

    let layout_method = matches
        .value_of("layout_method")
        .unwrap_or("2")
        .parse::<u8>()
        .unwrap();

    let bestn = matches
        .value_of("bestn")
        .unwrap_or("6")
        .parse::<usize>()
        .unwrap();

    let parameters = Parameters {
        nchunks: nchunks,
        nthreads: nthreads,
        w: wsize,
        k: ksize,
        r: rfactor,
        tol: tol,
        min_ec_cov: min_ec_cov,
    };

    log::info!("pg_asm {}", VERSION_STRING);
    log::info!(
        "command: {}",
        std::env::args().collect::<Vec<String>>().join(" ")
    );
    let cdir = std::env::current_dir()?;
    log::info!("current dir: {}", cdir.as_os_str().to_string_lossy());

    let start_wall_clock_time = SystemTime::now();
    log::info!("pg_asm run start");

    unsafe {
        log_resource("BGN: pg_asm", &mut rdata.assume_init_mut());
    }
    log::info!(
        "pg_asm run parameters: w:{}, k:{}, r:{}, tol:{} bestn:{}",
        wsize,
        ksize,
        rfactor,
        tol,
        bestn
    );

    log::info!("faster mode: {}", fastmode);
    log::info!("use layout method: {}", layout_method);
    log::info!("keep intermediate files: {}", keep);
    log::info!("number of threads: {}", nthreads);
    log::info!("number of chunks: {}", nchunks);
    log::info!("input read file: {}", input_reads);
    log::info!("working directory: {}", work_dir);
    log::info!(
        "sys: number of physical CPU cores detected: {}",
        physical_cpus
    );
    let system = sysinfo::System::new_all();
    log::info!("sys: total memory: {} KB", system.total_memory());
    log::info!("sys: used memory: {} KB", system.used_memory());
    log::info!("sys: total swap: {} KB", system.total_swap());
    log::info!("sys: used swap: {} KB", system.used_swap());

    if !Path::new(&work_dir).exists() {
        create_dir_all(&work_dir)?;
    };

    unsafe {
        get_ovlps(
            &input_reads,
            &work_dir,
            &prefix,
            false,
            &parameters,
            &mut rdata.assume_init_read(),
        )?;
    }

    if fastmode {
        log::info!("Fast mode: ignore read level error correction");
        // graph processing

        let seqdb = cat_path(&work_dir, &format!("{}.seqdb", &prefix));
        let seqidx = cat_path(&work_dir, &format!("{}.idx", &prefix));
        let ovlp_out = cat_path(&work_dir, &format!("{}-ovlp", &prefix));
        let layout_prefix = cat_path(&work_dir, &"asm".to_string());

        let layout_file = format!("{}_layout.dat", &layout_prefix);

        unsafe { log_resource("BGN: ovlp2layout", &mut rdata.assume_init_mut()); }
        log::info!("use layout method: {}", layout_method);
        match layout_method {
            1 => graph::ovlp2layout_v1(&ovlp_out, &layout_prefix, bestn),
            _ => dp_graph::ovlp2layout_v2(&ovlp_out, &layout_prefix, bestn)?,
        }
        unsafe { log_resource("END: ovlp2layout", &mut rdata.assume_init_mut()); }

        // layout -> sequence
        let output_file_prefix = format!("{}/asm_ctgs", &work_dir);
        unsafe { log_resource("BGN: layout2ctg", &mut rdata.assume_init_mut()); }
        layout::layout2ctg(&seqdb, &seqidx, &layout_file, &output_file_prefix)?;
        let _res = unsafe { getrusage(RUSAGE_SELF, &mut rdata.assume_init_read()) };
        unsafe { log_resource("END: layout2ctg", &mut rdata.assume_init_mut()); }
    } else {
        // error correction
        let seqdb = cat_path(&work_dir, &format!("{}.seqdb", &prefix));
        let seqidx = cat_path(&work_dir, &format!("{}.idx", &prefix));
        let shmmer_idx = cat_path(&work_dir, &format!("{}-shmr", &prefix));
        let ovlp_out = cat_path(&work_dir, &format!("{}-ovlp", &prefix));
        let ec_read_prefix = cat_path(&work_dir, &"ec_read".to_string());

        unsafe { log_resource("BGN: ovlp_ec", &mut rdata.assume_init_mut()); }
        ovlp_ec::ovlp_ec(&seqdb, &seqidx, &ovlp_out, &ec_read_prefix, &parameters)?;
        unsafe { log_resource("END: ovlp_ec", &mut rdata.assume_init_mut());}

        if !keep {
            let ovlp_out_ptn = format!("{}*", ovlp_out);
            for out in glob::glob(&ovlp_out_ptn.as_str()).expect("error to delete file") {
                let out = out.unwrap().to_string_lossy().into_owned();
                log::info!("remove {}", out);
                remove_file(out)?;
            }
            let shmmer_idx_ptn = format!("{}*", &shmmer_idx);
            for f in glob(&shmmer_idx_ptn.as_str()).expect("error to delete file") {
                let f = f.unwrap().to_string_lossy().into_owned();
                log::info!("remove {}", f);
                remove_file(f)?;
            }

            remove_file(&seqdb)?;
            log::info!("remove {}", seqdb.to_string());
            remove_file(&seqidx)?;
            log::info!("remove {}", seqidx.to_string());
        }

        let ec_lst = cat_path(&work_dir, &"ec_reads.lst".to_string());
        let mut ec_lst_file = BufWriter::new(File::create(&ec_lst).unwrap());
        let ec_file_ptn = format!("{}*.fa", &ec_read_prefix);
        for f in glob(&ec_file_ptn.as_str()).unwrap() {
            writeln!(ec_lst_file, "{}", f.unwrap().to_string_lossy())?;
        }
        drop(ec_lst_file); //close the file

        log::info!("input reads: {}", ec_lst);
        log::info!("working directory: {}", work_dir);

        let prefix = "ec_reads".to_string();
        if keep {
            unsafe { get_ovlps(&ec_lst, &work_dir, &prefix, false, &parameters, &mut rdata.assume_init_mut())?; }
        } else {
            unsafe { get_ovlps(&ec_lst, &work_dir, &prefix, true, &parameters, &mut rdata.assume_init_mut())?; }
        }

        if !keep {
            let ec_file_ptn = format!("{}*.fa", &ec_read_prefix);
            for f in glob(&ec_file_ptn.as_str()).unwrap() {
                let f = f.unwrap().to_string_lossy().into_owned();
                log::info!("remove {}", f);
                remove_file(f)?;
            }
            let shmmer_idx = cat_path(&work_dir, &format!("{}-shmr", &prefix));
            let shmmer_idx_ptn = format!("{}*", &shmmer_idx);
            for f in glob(&shmmer_idx_ptn.as_str()).expect("error to delete file") {
                let f = f.unwrap().to_string_lossy().into_owned();
                log::info!("remove {}", f);
                remove_file(f)?;
            }
        }

        //step 4: graph processing

        let seqdb = cat_path(&work_dir, &format!("{}.seqdb", &prefix));
        let seqidx = cat_path(&work_dir, &format!("{}.idx", &prefix));
        let ovlp_out = cat_path(&work_dir, &format!("{}-ovlp", &prefix));
        let layout_prefix = cat_path(&work_dir, &"asm".to_string());
        let layout_file = format!("{}_layout.dat", &layout_prefix);

        unsafe {log_resource("BGN: ovlp2layout", &mut rdata.assume_init_mut());}
        log::info!("use layout method: {}", layout_method);
        match layout_method {
            1 => graph::ovlp2layout_v1(&ovlp_out, &layout_prefix, bestn),
            _ => dp_graph::ovlp2layout_v2(&ovlp_out, &layout_prefix, bestn)?,
        }
        unsafe{log_resource("END: ovlp2layout", &mut rdata.assume_init_mut());}

        //step 5: layout -> sequence
        unsafe{log_resource("BGN: layout2ctg", &mut rdata.assume_init_mut());}
        let output_file_prefix = format!("{}/asm_ctgs", &work_dir);
        layout::layout2ctg(&seqdb, &seqidx, &layout_file, &output_file_prefix)?;
        let _res = unsafe { getrusage(RUSAGE_SELF, &mut rdata.assume_init_read()) };
        unsafe{log_resource("END: layout2ctg", &mut rdata.assume_init_mut());}
    }
    if no_resolve {
        log::info!("ignore dup resolution");
    } else {
        let ref_file = format!("{}/asm_ctgs_m.fa", &work_dir);
        let tgt_file = format!("{}/asm_ctgs_e0.fa", &work_dir);
        let out_file = format!("{}/asm_ctgs_e.fa", &work_dir);

        unsafe{log_resource("BEN: dedup_a_ctgs", &mut rdata.assume_init_mut());}
        dedup_target_seqs(&ref_file, &tgt_file, &out_file, wsize, ksize, rfactor)?;
        unsafe{log_resource("END: dedup_a_ctgs", &mut rdata.assume_init_mut());}

        let resolve_prefix = format!("{}/asm_ctgs_m", &work_dir);

        unsafe{log_resource("BEN: resolve_ht", &mut rdata.assume_init_mut());}
        resolve_ht(&ref_file, &resolve_prefix, wsize, ksize, rfactor)?;
        unsafe{log_resource("END: resolve_ht", &mut rdata.assume_init_mut());}
    }
    let (_, ut, st) = unsafe { log_resource("END: pg_asm", &mut rdata.assume_init_mut()) };
    log::info!("pg_asm run end");
    log::info!(
        "total user cpu time: {} seconds = {} hours",
        ut,
        ut as f32 / 60.0 / 60.0
    );
    log::info!(
        "total system cpu time: {} seconds = {} hours",
        st,
        st as f32 / 60.0 / 60.0
    );
    let elapsed_time = start_wall_clock_time.elapsed().unwrap().as_secs_f32();
    log::info!(
        "total elapse time: {} seconds = {} hours",
        elapsed_time,
        elapsed_time / 60.0 / 60.0
    );
    Ok(())
}
