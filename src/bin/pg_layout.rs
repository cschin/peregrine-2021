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
use utils::layout::layout2ctg;

fn main() -> Result<(), std::io::Error> {
    let matches = clap_app!(pg_layout =>
        (version: VERSION_STRING)
        (author: "Jason Chin <jason@omnibio.ai>")
        (about: "
Peregrine-2021 genome assembler, 
pg_layout: convert the assembly graph to paths and generate the contig fasta file
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
        (@arg SEQDB:+required "Path to the seqdb file ")
        (@arg SEQIDX:+required "Path to the seqdb index file")
        (@arg layout_file: --layout_file +required +takes_value "Path to the layout file")
        (@arg output_file: --output +required +takes_value "Path to the output file")
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
    let layout_file = matches.value_of("layout_file").unwrap().to_string();
    let output_file = matches.value_of("output_file").unwrap().to_string();

    let _res = log::info!("layout:seq_db_file: {}", seqdb_file);
    let _res = log::info!("layout:index_file: {}", index_file);
    let _res = log::info!("layout:layout_file: {}", layout_file);
    let _res = log::info!("layout:output: {}", output_file);

    layout2ctg(&seqdb_file, &index_file, &layout_file, &output_file)?;
    Ok(())
}
