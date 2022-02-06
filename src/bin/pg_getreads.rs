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
use memmap::MmapOptions;
use simple_logger::SimpleLogger;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
mod utils;
use utils::shmmrutils::{get_seq_fragment, ReadLocation};

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn main() -> () {
    let matches = clap_app!(pg_getreads =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler, 
pg_getreads: generate fasta file for a subset of reads from the sequence database
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg SEQDB: +required "Path to the seqdb file ")
        (@arg SEQIDX: +required "Path to the seqdb index file")
        (@arg READID: +required "Path to the read id file")
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

    let seqdb_file = matches.value_of("SEQDB").unwrap().to_string();
    let index_file = matches.value_of("SEQIDX").unwrap().to_string();
    let read_id_file = matches.value_of("READID").unwrap().to_string();

    let _res = writeln!(io::stderr(), "seq_db_file: {}", seqdb_file);
    let _res = writeln!(io::stderr(), "index_file: {}", index_file);
    let _res = writeln!(io::stderr(), "read_id_file: {}", read_id_file);

    let mut read_index = Vec::<ReadLocation>::new();
    let mut read_name = Vec::<String>::new();

    if let Ok(lines) = read_lines(index_file) {
        for line in lines {
            if let Ok(rec) = line {
                //let rec_trimmed = rec.trim_end();
                // the record line looks like 000000023 m64062_190803_042216/144/ccs 20359 467415
                let v: Vec<&str> = rec.split_whitespace().collect();
                //let rid: u32 = v[0].parse().unwrap();
                let start: usize = v[3].parse().unwrap();
                let len: usize = v[2].parse().unwrap();
                read_index.push(ReadLocation {
                    start: start,
                    len: len,
                });
                read_name.push(v[1].to_string());
                //println!("{} {} {}", rid, start, len);
            }
        }
    }

    let file = File::open(seqdb_file).unwrap();
    let mmap = unsafe { MmapOptions::new().map(&file).unwrap() };
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    if let Ok(lines) = read_lines(read_id_file) {
        for line in lines {
            if let Ok(rec) = line {
                let rec_trimmed = rec.trim_end();
                let rid0 = rec_trimmed.parse::<u32>().unwrap();
                let rloc = read_index[rid0 as usize];
                let len = rloc.len as u32;
                let seq_frag = get_seq_fragment(rid0, 0, 0, len, &mmap, &read_index);
                let _ = writeln!(handle, ">{} {:09}", read_name[rid0 as usize], rid0);
                let _ = writeln!(handle, "{}", String::from_utf8_lossy(&seq_frag));
            }
        }
    }
}
