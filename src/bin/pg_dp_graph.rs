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
use utils::dp_graph::ovlp2layout_v2;
fn main() -> Result<(), std::io::Error> {
    let matches = clap_app!(pg_graph =>
            (version: VERSION_STRING)
            (author: "Jason Chin <jason@omnibio.ai>")
            (about: "
Peregrine-2021 genome assembler, 
pg_dp_graph: take overlap data file as input to generate the layout file using a polyploid aware layout algorithm
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/")    
            (@arg prefix: --prefix +required +takes_value "Path prefix for input files")
            (@arg out_prefix: --out_prefix +required +takes_value "ath prefix for output files ")
            (@arg bestn: --bestn -b +takes_value "number of best overlaps for initial graph [default: 6]")
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
        
    let prefix = matches.value_of("prefix").unwrap().to_string();
    let out_prefix = matches.value_of("out_prefix").unwrap().to_string();
    let bestn = matches
        .value_of("bestn")
        .unwrap_or("6")
        .parse::<usize>()
        .unwrap();

    let _err = log::info!("graph:out_prefix: {}", out_prefix,);
    ovlp2layout_v2(&prefix, &out_prefix, bestn)?;
    Ok(())
}
