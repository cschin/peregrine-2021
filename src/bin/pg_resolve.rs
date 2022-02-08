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
use utils::resolve::resolve_ht;

fn main() -> Result<(), std::io::Error> {
    let matches = clap_app!(pg_resolve =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler, 
pg_resolve: this tool aligns all contigs to themselve to identify haplotype-related contigs
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg fasta_file: -f --fasta_file +required +takes_value "Path to the layout file")
        (@arg output_prefix: -o --out_prefix +required +takes_value "Path to the output prefix")
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

    let fasta_file = matches.value_of("fasta_file").unwrap().to_string();
    let output_prefix = matches.value_of("output_prefix").unwrap().to_string();
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

    resolve_ht(&fasta_file, &output_prefix, wsize, ksize, rfactor)?;
    Ok(())
}
