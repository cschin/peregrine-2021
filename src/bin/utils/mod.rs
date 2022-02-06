// Peregrine Assembler and SHIMMER Genome Assembly Toolkit 
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the 
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

pub mod build_idx;
pub mod build_sdb;
pub mod dp_graph;
pub mod seqmap;
pub mod graph;
pub mod graph_analysis;
pub mod layout;
pub mod ovlp;
pub mod ovlp_ec;
pub mod resolve;
pub mod shmmrutils;
pub use core::mem::MaybeUninit;
pub use libc::{getrusage, rusage, RUSAGE_SELF, RUSAGE_THREAD};

#[derive(Copy, Clone)]
pub struct Parameters {
    pub nthreads: u32,
    pub nchunks: u32,
    pub k: u32,
    pub w: u32,
    pub r: u32,
    pub tol: f64,
    pub min_ec_cov: u16,
}

#[allow(dead_code)]
pub fn log_resource(msg: &str, data: &mut rusage) -> (u64, u64, u64) {
    let _res = unsafe { getrusage(RUSAGE_SELF, data) };
    log::info!(
        "{} : (maxRSS, utime, stime): {} {} {}",
        msg,
        data.ru_maxrss,
        data.ru_utime.tv_sec,
        data.ru_stime.tv_sec
    );

    (
        data.ru_maxrss as u64,
        data.ru_utime.tv_sec as u64,
        data.ru_stime.tv_sec as u64,
    )
}
