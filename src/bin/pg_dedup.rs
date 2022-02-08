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
use simple_logger::SimpleLogger;

mod utils;
use utils::seqmap;

fn main() -> Result<(), std::io::Error> {
    let matches = clap_app!(pg_resolve =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler, 
pg_dedup: perform all contigs to all contigs alignment to remove duplicates 
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg ref_fasta: -f --ref_fasta +required +takes_value "Path to the reference file")
        (@arg target_fasta: -t --target_fasta +required +takes_value "Path to the target file")
        (@arg output: -o --output +required +takes_value "Path to the output filename")
        (@arg w: -w +takes_value "Window size [default: 48]")
        (@arg k: -k +takes_value "Kmer size [default: 56]")
        (@arg r: -r +takes_value "Reduction factor [default: 4]")
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

    let ref_fasta_file = matches.value_of("ref_fasta").unwrap().to_string();
    let target_fasta_file = matches.value_of("target_fasta").unwrap().to_string();
    let output_file = matches.value_of("output").unwrap().to_string();
    let wsize = matches
        .value_of("w")
        .unwrap_or("48")
        .parse::<u32>()
        .unwrap();

    let ksize = matches
        .value_of("k")
        .unwrap_or("56")
        .parse::<u32>()
        .unwrap();

    let rfactor = matches.value_of("r").unwrap_or("4").parse::<u32>().unwrap();
    seqmap::dedup_target_seqs(
        &ref_fasta_file,
        &target_fasta_file,
        &output_file,
        wsize,
        ksize,
        rfactor,
    )?;
    Ok(())
}
