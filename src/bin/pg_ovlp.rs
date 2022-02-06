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
mod utils;
use simple_logger::SimpleLogger;
use utils::ovlp::ovlp;
use utils::Parameters;
fn main() -> Result<(), std::io::Error> {
    let matches = clap_app!(pg_ovlp =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (@arg SEQDB:+required "Path to the seqdb file ")
        (about: "
Peregrine-2021 genome assembler, 
pg_ovlp: generate haplotype specific overlaps between the reads
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg SEQIDX:+required "Path to the seqdb index file")
        (@arg SHMRINDEXPREFIX: +required "The prefix to the output shimmer index database")
        (@arg OUTPREFIX: +required "The prefix of the output ovelap files")
        (@arg NTHREADS: +required "Number of threads ")
        (@arg NCHUNKS: +required "Number of partition")
        (@arg tol: -t --tol +takes_value "Alignment tolerance [default: 0.01]")
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
    let shimmer_index_file_prefix = matches.value_of("SHMRINDEXPREFIX").unwrap().to_string();
    let out_prefix = matches.value_of("OUTPREFIX").unwrap().to_string();
    let nthreads = matches
        .value_of("NTHREADS")
        .unwrap()
        .parse::<u32>()
        .unwrap();
    let nchunks = matches.value_of("NCHUNKS").unwrap().parse::<u32>().unwrap();
    let tol = matches
        .value_of("tol")
        .unwrap_or("0.01")
        .parse::<f64>()
        .unwrap();

    let parameters = Parameters {
        nchunks: nchunks,
        nthreads: nthreads,
        w: 0,
        k: 0,
        r: 0,
        tol: tol,
        min_ec_cov: 1,
    };

    ovlp(
        &seqdb_file,
        &index_file,
        &shimmer_index_file_prefix,
        &out_prefix,
        &parameters,
    )?;
    Ok(())
}
