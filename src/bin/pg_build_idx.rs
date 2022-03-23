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
use utils::build_idx::build;
use utils::Parameters;
fn main() -> () {
    let matches = clap_app!(pg_build_idx =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler
build the SHIMMER index from the reads for overlapping
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")            
        (@arg SEQDB:+required "Path to the seqdb file ")
        (@arg SEQIDX:+required "Path to the seqdb index file")
        (@arg SHMRINDEXPREFIX: +required "The prefix to the output shimmer index database")
        (@arg NTHREADS: +required "Number of threads")
        (@arg NCHUNKS: +required "Number of partition")
        (@arg w: -w +takes_value "Window size [default: 80]")
        (@arg k: -k +takes_value "Kmer size [default: 56]")
        (@arg r: -r +takes_value "Reduction factor [default: 6]")
        (@arg log: --log +takes_value "log level: DBBUG or INFO (default)")
    )
    .get_matches();

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
    let out_prefix = matches.value_of("SHMRINDEXPREFIX").unwrap().to_string();
    let nthreads = matches
        .value_of("NTHREADS")
        .unwrap()
        .parse::<u32>()
        .unwrap();
    let nchunks = matches.value_of("NCHUNKS").unwrap().parse::<u32>().unwrap();
    let wsize = matches.value_of("w").unwrap_or("80").parse::<u32>().unwrap();
    let ksize = matches.value_of("k").unwrap_or("56").parse::<u32>().unwrap();
    let rfactor = matches.value_of("r").unwrap_or("6").parse::<u32>().unwrap();

    let parameters = Parameters {
        nchunks: nchunks,
        nthreads: nthreads,
        w: wsize,
        k: ksize,
        r: rfactor,
        tol: 0.0, //not used
        min_ec_cov: 1,
    };

    build(&seqdb_file, &index_file, &out_prefix, &parameters);
}
