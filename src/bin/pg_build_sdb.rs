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
use utils::build_sdb::build;
fn main() -> () {
    let matches = clap_app!(pg_build_sdb =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler 
pb_build_sdb: build the sequence database
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg input: --input +required +takes_value "Path to a file that contains the list of reads in .fa .fa.gz .fastq or fastq.gz formats")
        (@arg out_prefix: --out_prefix +required +takes_value "The prefix for the sequence database and index files")
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

    let seq_list_file = matches.value_of("input").unwrap().to_string();
    let out_prefix = matches.value_of("out_prefix").unwrap().to_string();
    let _nbases = build(&seq_list_file, &out_prefix);
}
