// Peregrine Assembler and SHIMMER Genome Assembly Toolkit
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

#![allow(dead_code)]

//
// utility functions for SHIMMER process and the read matching / mapping
//

use memmap::Mmap;
use std::fmt;

use std::collections::HashMap;

pub struct HPCSeq {
    pub s: Vec<u8>,
    pub p: Vec<u32>,
}

pub fn get_hpc_seq(seq0: &Vec<u8>) -> HPCSeq {
    let mut seq = HPCSeq {
        s: Vec::<u8>::with_capacity(seq0.len()),
        p: Vec::<u32>::with_capacity(seq0.len()),
    };
    let mut i = 0_usize;
    let seq0len = seq0.len();

    while i < seq0len - 1 {
        let b = seq0[i]; // 2bit flag + 2bit base
        if b & 0b1100 != 0 {
            i += 1; // base marked as hp/dimer tial
            continue;
        }
        seq.s.push(b & 0b0011);
        seq.p.push(i as u32);
        i += 1;
    }
    seq
}

pub struct OvlpMatch {
    pub m_size: u32,
    pub dist: u32,
    pub bgn0: u32,
    pub end0: u32,
    pub bgn1: u32,
    pub end1: u32,
    pub m_end0: u32,
    pub m_end1: u32,
    pub deltas: Option<Vec<DeltaPoint>>,
}

#[derive(Copy, Clone)]
pub struct DeltaPoint {
    pub x: u32,
    pub y: u32,
    pub dk: i32,
}

pub fn track_delta_point(
    delta_pts: &HashMap<(u32, i32), DeltaPoint>,
    d_final: u32,
    k_final: i32,
    s: u32,
    e: u32,
) -> Vec<DeltaPoint> {
    let mut dpts = Vec::<DeltaPoint>::with_capacity(d_final as usize);
    let mut d = d_final;
    let mut k = k_final;
    while d > 0 {
        let dpt = delta_pts.get(&(d, k)).unwrap();
        if dpt.x > s && dpt.x < e {
            dpts.push(*dpt);
        }
        d -= 1;
        k -= dpt.dk;
    }
    dpts
}

pub fn match_reads(
    seq0: &Vec<u8>,
    seq1: &Vec<u8>,
    get_delta: bool,
    tol: f64,
    min_match_len: u32,
    bandwidth: u32,
) -> Option<OvlpMatch> {
    //
    // A variation of the O(nD) algorithm for read alignments
    //

    // let min_match_len = 1200;
    let len0 = seq0.len();
    let len1 = seq1.len();
    //println!("S {} {}", len0, len1);
    //let d_max = 64 + (0.01 * if len0 < len1 {len0 as f32} else {len1 as f32}) as u32;
    let d_max = 32
        + (tol
            * if len0 < len1 {
                len0 as f64
            } else {
                len1 as f64
            }) as u32;
    let max_band_width = bandwidth;
    let band_tolerance = bandwidth;
    let mut k_min = 0_i32;
    let mut k_max = 0_i32;
    let mut uv_map = HashMap::<i32, (u32, u32)>::new();
    // uv_map: maping k to the u, v, which keep the d path end in k
    let mut delta_pts = HashMap::<(u32, i32), DeltaPoint>::new();
    let mut x: u32;
    let mut y: u32;
    let mut x1: u32;
    let mut y1: u32;
    let mut best_m = -1_i32;
    let mut matched = false;
    let mut d_final = 0_u32;
    let mut k_final = 0_i32;
    let mut pre_k: i32;
    let mut start = false;
    let mut longest_match = 0_u32;
    let mut rtn = OvlpMatch {
        m_size: 0,
        dist: 0,
        bgn0: 0,
        end0: 0,
        bgn1: 0,
        end1: 0,
        m_end0: 0,
        m_end1: 0,
        deltas: None,
    };

    for d in -(d_max as i32)..=(d_max as i32) {
        uv_map.insert(d, (0, 0));
    }
    for d in 0..d_max {
        if k_max - k_min > max_band_width as i32 {
            // println!("KK {} {} {} {}", k_max, k_min, k_max - k_min, max_band_width);
            break;
        }
        for k in (k_min..=k_max).step_by(2) {
            let (_, vn) = uv_map.get(&(k - 1)).unwrap();
            let (_, vp) = uv_map.get(&(k + 1)).unwrap();
            if k == k_min || ((k != k_max) && vn < vp) {
                x = *vp;
                pre_k = k + 1;
            } else {
                x = *vn + 1;
                pre_k = k - 1;
            }
            y = ((x as i32) - k) as u32;

            if get_delta {
                let dpt = DeltaPoint {
                    x: x,
                    y: y,
                    dk: k - pre_k,
                };
                delta_pts.entry((d, k)).or_insert(dpt);
            };

            x1 = x;
            y1 = y;

            while (x as usize) < len0 - 1
                && (y as usize) < len1 - 1
                && seq0[x as usize] == seq1[y as usize]
            {
                x += 1;
                y += 1;
            }

            if (x - x1) > 8 {
                if !start {
                    rtn.bgn0 = x1;
                    rtn.bgn1 = y1;
                    start = true;
                }
                // we set the ends here to avoid bad sequences
                // this way, we are sure that, at least, 8 bases are aligned
                rtn.end0 = x;
                rtn.end1 = y;
            }

            if (x - x1) > longest_match {
                longest_match = x - x1;
                rtn.m_end0 = x;
                rtn.m_end1 = y;
            }

            // println!("IM {} {} {} {} {} {} {} {}", x, y, len0, len1, d, d_max, k, pre_k);
            uv_map.insert(k, (x + y, x));
            if (x + y) as i32 > best_m {
                best_m = (x + y) as i32;
            }
            if (x as usize) >= len0 - 1 || (y as usize) >= len1 - 1 {
                matched = true;
                d_final = d;
                k_final = k;
                break;
            }
        }
        // For banding
        let mut k_max_new = k_min;
        let mut k_min_new = k_max;
        for k2 in (k_min..=k_max).step_by(2) {
            let (u, _) = uv_map.get(&k2).unwrap();
            if *u as i32 >= (best_m - (band_tolerance as i32)) {
                if k2 < k_min_new {
                    k_min_new = k2;
                }
                if k2 > k_max_new {
                    k_max_new = k2;
                }
            }
        }

        k_max = k_max_new + 1;
        k_min = k_min_new - 1;
        if matched == true {
            //println!("match: {} {}", d_final, k_final);
            let mut d_inside = 0_u32;
            if get_delta {
                let dpts = track_delta_point(&delta_pts, d_final, k_final, rtn.bgn0, rtn.end0);
                for dpt in &dpts {
                    if dpt.x > rtn.bgn0 && dpt.x < rtn.end0 {
                        d_inside += 1;
                    }
                }
                rtn.deltas = Some(dpts);
            }
            rtn.dist = d_inside;
            rtn.m_size = (rtn.end0 - rtn.bgn0 + rtn.end1 - rtn.bgn1 + 2 * d_inside) >> 1;
            if rtn.m_size < min_match_len {
                matched = false;
            }
            break;
        }
    }
    if !matched {
        None
    } else {
        Some(rtn)
    }
}

#[derive(Clone, Copy)]
pub struct ReadLocation {
    pub start: usize,
    pub len: usize,
}

pub fn get_2bit_fragment(
    rid: u32,
    strand: u8,
    bgn: u32,
    end: u32,
    readsdb: &Mmap,
    read_index: &Vec<ReadLocation>,
) -> Vec<u8> {
    let mut seq = Vec::<u8>::new();
    let rloc = read_index[rid as usize];
    let s = rloc.start + bgn as usize;
    let e = rloc.start + end as usize;

    if strand == 0 {
        for c in &readsdb[s..e] {
            seq.push(c & 0b0011);
        }
    } else {
        for c in &readsdb[s..e] {
            seq.push((c >> 4) & 0b0011);
        }
    }
    seq
}

pub fn get_seq_fragment(
    rid: u32,
    strand: u8,
    bgn: u32,
    end: u32,
    mmap: &Mmap,
    read_index: &Vec<ReadLocation>,
) -> Vec<u8> {
    let mut seq = Vec::<u8>::new();
    let base_map = &[b'A', b'C', b'G', b'T'];
    let rloc = read_index[rid as usize];
    let s = rloc.start + bgn as usize;
    let e = rloc.start + end as usize;

    if strand == 0 {
        for c in &mmap[s..e] {
            if c & 0b1100 != 0b1100 {
                seq.push(base_map[(c & 0b0011) as usize]);
            } else {
                seq.push(b'N');
            }
        }
    } else {
        for c in &mmap[s..e] {
            if ((c >> 4) & 0b1100) != 0b1100 {
                seq.push(base_map[((c >> 4) & 0b0011) as usize]);
            } else {
                seq.push(b'N');
            }
        }
    }
    //println!("{} {}", rid, String::from_utf8_lossy(&seq));
    seq
}

#[derive(Clone, Copy)]
pub struct MM128 {
    pub x: u64,
    pub y: u64,
}

impl fmt::Display for MM128 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let hash = self.x >> 8;
        let span = (self.x & 0xFF) as u8;
        let rid = (self.y >> 32) as u32;
        let pos = ((self.y & 0xFFFFFFFF) >> 1) as u32;
        let strand = (self.y & 0x1) as u8;
        write!(f, "({}, {}, {}, {}, {})", hash, span, rid, pos, strand)
    }
}

pub fn u64hash(key: u64) -> u64 {
    let key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    let key = key ^ key >> 24;
    let key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    let key = key ^ key >> 14;
    let key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    let key = key ^ key >> 28;
    let key = key.wrapping_add(key << 31);
    key
}

fn _u64hash(key: u64) -> u64 {
    let key = !key + (key << 21); // key = (key << 21) - key - 1;
    let key = key ^ key >> 24;
    let key = (key + (key << 3)) + (key << 8); // key * 265
    let key = key ^ key >> 14;
    let key = (key + (key << 2)) + (key << 4); // key * 21
    let key = key ^ key >> 28;
    let key = key + (key << 31);
    key
}

pub struct RingBuffer {
    pub v: Vec<MM128>,
    pub size: usize,
    pub start_pos: usize,
    pub end_pos: usize,
    pub len: usize,
}

impl RingBuffer {
    pub fn new(size: usize) -> Self {
        let vv = vec![
            MM128 {
                x: u64::MAX,
                y: u64::MAX
            };
            size
        ];
        RingBuffer {
            v: vv,
            size: size,
            start_pos: 0,
            end_pos: 0,
            len: 0,
        }
    }

    pub fn push(&mut self, m: MM128) {
        if self.len < self.size {
            self.v[self.end_pos] = m;
            self.end_pos += 1;
            self.end_pos %= self.size;
            self.len += 1;
        } else {
            self.v[self.end_pos] = m;
            self.end_pos += 1;
            self.end_pos %= self.size;
            self.start_pos += 1;
            self.start_pos %= self.size;
        }
    }

    pub fn _clear(&mut self) {
        self.v.clear();
        self.len = 0;
        self.start_pos = 0;
        self.end_pos = 0;
    }

    pub fn get_min(&self) -> MM128 {
        let mut min = MM128 {
            x: u64::MAX,
            y: u64::MAX,
        };
        for i in 0..self.len {
            if self.v[i].x < min.x {
                min = self.v[i];
            }
        }
        min
    }

    pub fn get(&self, i: usize) -> MM128 {
        self.v[(self.start_pos + i) % self.size]
    }
}

pub fn reduce_shmmr(mers: Vec<MM128>, r: u32) -> Vec<MM128> {
    let mut shmmrs = Vec::<MM128>::new();
    let mut rbuf = RingBuffer::new(r as usize);
    let mut min_mer = MM128 {
        x: u64::MAX,
        y: u64::MAX,
    };
    let mut pos = 0;
    let mut mdist = 0;
    loop {
        if pos >= mers.len() {
            break;
        }
        let m = mers[pos];
        rbuf.push(m);
        if mdist == (r - 1) as usize {
            min_mer = rbuf.get_min();
            let mut last_i = 0_usize;
            for i in 0..rbuf.size as usize {
                let mm = rbuf.get(i);
                if mm.x == min_mer.x {
                    shmmrs.push(mm);
                    min_mer = mm;
                    last_i = i;
                }
            }
            mdist = r as usize - 1 - last_i;
            pos += 1;
            continue;
        } else if m.x <= min_mer.x && pos >= r as usize {
            shmmrs.push(m);
            min_mer = m;
            mdist = 0;
            pos += 1;
            continue;
        }
        mdist += 1;
        pos += 1;
    }
    shmmrs
}

pub fn sequence_to_shmmrs(rid: u32, seq: &Vec<u8>, w: u32, k: u32, r: u32) -> Vec<MM128> {
    //let base2bits: [u64; 4] = [0, 1, 2, 3];

    let mut shmmrs = Vec::<MM128>::new();

    let mut pos = 0;
    let mut mdist = 0;
    let shift = k - 1;
    assert!(k <= 56);
    assert!(w <= 128);
    assert!(r > 0 && r < 13);
    let mut fmmer = (0_u64, 0_u64);
    let mut rmmer = (0_u64, 0_u64);
    let mask = u64::MAX >> (64 - k);
    let mut rbuf = RingBuffer::new(w as usize);
    let mut min_mer = MM128 {
        x: u64::MAX,
        y: u64::MAX,
    };
    loop {
        if pos >= seq.len() {
            break;
        }
        let c = (seq[pos] & 0b0011) as u64;
        // println!("C {} {} {}", seq[pos], pos, c);
        if seq[pos] & 0b1100 != 0b1100 {
            // Not non-A,C,G,T base
            fmmer.0 <<= 1;
            fmmer.0 |= c & 0b01;
            fmmer.0 &= mask;
            fmmer.1 <<= 1;
            fmmer.1 |= (c & 0b10) >> 1;
            fmmer.1 &= mask;

            let rc = 0x03 ^ c;
            rmmer.0 >>= 1;
            rmmer.0 |= (rc & 0b01) << shift;
            rmmer.0 &= mask;
            rmmer.1 >>= 1;
            rmmer.1 |= ((rc & 0b10) >> 1) << shift;
            rmmer.1 &= mask;
        }
        if fmmer == rmmer {
            pos += 1;
            continue;
        }
        if pos < k as usize {
            pos += 1;
            continue;
        }
        let mut forward = true;
        if rmmer.0 < fmmer.0 {
            forward = false;
        }
        let mmer_hash = match forward {
            true => u64hash(fmmer.0) ^ u64hash(fmmer.1) ^ 0x0,
            false => u64hash(rmmer.0) ^ u64hash(rmmer.1) ^ 0x0,
        };
        let strand: u64 = if forward { 0 } else { 1 };
        let m = MM128 {
            x: mmer_hash << 8 | k as u64,
            y: (rid as u64) << 32 | (pos as u64) << 1 | strand,
        };
        rbuf.push(m);
        // println!("MM {} {} {} {}", rid, m, fmmer.0, fmmer.1);
        if mdist == (w - 1) as usize {
            min_mer = rbuf.get_min();
            for i in 0..rbuf.size as usize {
                let mm = rbuf.get(i);
                if mm.x == min_mer.x {
                    shmmrs.push(mm);
                    min_mer = mm;
                    // println!("MM0 {} {}", rid, mm);
                }
            }
            mdist = pos - ((min_mer.y & 0xFFFFFFFF) >> 1) as usize;
            pos += 1;
            continue;
        } else if m.x <= min_mer.x && pos >= w as usize {
            shmmrs.push(m);
            // println!("MM1 {} {}", rid, m);
            min_mer = m;
            mdist = 0;
            pos += 1;
            continue;
        }
        mdist += 1;
        pos += 1;
    }
    let shmmrs = reduce_shmmr(shmmrs, r);
    let shmmrs = reduce_shmmr(shmmrs, r);
    shmmrs
}
