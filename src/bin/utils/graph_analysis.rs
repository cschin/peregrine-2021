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
// define the overlap and graph data structure and lower level graph processing utility functions
//

use petgraph::graphmap::DiGraphMap;
use petgraph::visit::Bfs;
use petgraph::Direction::{Incoming, Outgoing};
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;

use std::fs::File;
use std::io::{self, BufWriter, Write};

pub type U32AsmGraph = DiGraphMap<(u32, u8), u32>;
pub type OvlpGraph = U32AsmGraph;
pub type UtgGraph = U32AsmGraph;
pub type ReadNode = (u32, u8);
pub type OvlpEdge = (ReadNode, ReadNode);

#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
pub struct ReadPair {
    pub rid0: u32,
    pub strand0: u8,
    pub rid1: u32,
    pub strand1: u8,
}

impl ReadPair {
    pub fn new(v: (u32, u8), w: (u32, u8)) -> Self {
        ReadPair {
            rid0: v.0,
            strand0: v.1,
            rid1: w.0,
            strand1: w.1,
        }
    }
    pub fn _to_str(&self) -> String {
        format!(
            "{} {} {} {}",
            self.rid0, self.strand0, self.rid1, self.strand1
        )
    }

    pub fn reverse(&self) -> ReadPair {
        ReadPair {
            rid0: self.rid1,
            strand0: 1 - self.strand1,
            rid1: self.rid0,
            strand1: 1 - self.strand0,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Overlap {
    // main struct for keeping read overlap
    pub rid0: u32,
    pub rid1: u32,
    pub strand1: u8,
    pub len0: u32,
    pub len1: u32,
    pub d_left: i32,
    pub d_right: i32,
    pub bgn0: u32,
    pub end0: u32,
    pub bgn1: u32,
    pub end1: u32,
    pub dist: u32,
     // dist: raw distice determine by the O(Nd) alignment algorithm
    pub idt: f32,
    pub dist_c: u32,
    // dist_c: "distance" after hp corrections
    pub max_dist_c: u32,
    pub idt_c: f32,
    pub flag: u8,
    // flag bit field, Not used now 2020/11/02
    // 0x01: the rid0 is chimer
    // 0x02: the rir0 and rid1 are compatitable pair
    // 0x04: the rid1 is the best right pair
    // 0x08: the rid1 is the best left pair
    // 0x10: the rid1 is a chimer
    // 0x20: the rid1 is contained
    // 0x40: the rid0 is contained
}

impl Overlap {
    pub fn _new() -> Self {
        Self {
            rid0: 0,
            rid1: 0,
            strand1: 0,
            len0: 0,
            len1: 0,
            d_left: 0,
            d_right: 0,
            bgn0: 0,
            end0: 0,
            bgn1: 0,
            end1: 0,
            dist: 0,
            idt: 0.0,
            dist_c: 0,
            max_dist_c: 0,
            idt_c: 0.0,
            flag: 0,
        }
    }

    pub fn build_from(v: Vec<&str>) -> Self {
        Overlap {
            rid0: v[1].parse().unwrap(),
            rid1: v[2].parse().unwrap(),
            strand1: v[3].parse().unwrap(),
            len0: v[4].parse().unwrap(),
            len1: v[5].parse().unwrap(),
            d_left: v[6].parse().unwrap(),
            d_right: v[7].parse().unwrap(),
            bgn0: v[8].parse().unwrap(),
            end0: v[9].parse().unwrap(),
            bgn1: v[10].parse().unwrap(),
            end1: v[11].parse().unwrap(),
            dist: v[12].parse().unwrap(),
            idt: v[13].parse().unwrap(),
            dist_c: v[14].parse().unwrap(),
            max_dist_c: v[15].parse().unwrap(),
            idt_c: v[16].parse().unwrap(),
            flag: v[17].parse().unwrap(),
        }
    }

    pub fn _format(&self) -> String {
        format!(
            "{} {} {} {} {} {} {} {} {} {} {} {} {:.2} {} {} {:.2} {}",
            self.rid0,
            self.rid1,
            self.strand1,
            self.len0,
            self.len1,
            self.d_left,
            self.d_right,
            self.bgn0,
            self.end0,
            self.bgn1,
            self.end1,
            self.dist,
            self.idt,
            self.dist_c,
            self.max_dist_c,
            self.idt_c,
            self.flag
        )
    }

    pub fn swap_rp(&self) -> Overlap {
        // swap the overlapped pair
        let d_left: i32;
        let d_right: i32;
        let strand1: u8;
        let bgn0: u32;
        let end0: u32;
        let bgn1: u32;
        let end1: u32;
        if self.strand1 == 0 {
            d_left = -self.d_left;
            d_right = -self.d_right;
            strand1 = 0;
            bgn0 = self.bgn0;
            end0 = self.end0;
            bgn1 = self.bgn1;
            end1 = self.end1;
        } else {
            d_left = self.d_right;
            d_right = self.d_left;
            strand1 = 1;
            bgn0 = self.len0 - self.end0;
            end0 = self.len0 - self.bgn0;
            bgn1 = self.len1 - self.end1;
            end1 = self.len1 - self.bgn1;
        }
        Overlap {
            rid0: self.rid1,
            rid1: self.rid0,
            strand1: strand1,
            len0: self.len1,
            len1: self.len0,
            d_left: d_left,
            d_right: d_right,
            bgn0: bgn1,
            end0: end1,
            bgn1: bgn0,
            end1: end0,
            dist: self.dist,
            idt: self.idt,
            dist_c: self.dist_c,
            max_dist_c: self.max_dist_c,
            idt_c: self.idt_c,
            flag: self.flag,
        }
    }

    pub fn reverse_strand(&self) -> Overlap {
        // reverse the overlapped strain of the read 1
        Overlap {
            rid0: self.rid0,
            rid1: self.rid1,
            strand1: 1 - self.strand1,
            len0: self.len1,
            len1: self.len0,
            d_left: -self.d_right,
            d_right: -self.d_left,
            bgn0: self.len0 - self.end0,
            end0: self.len0 - self.bgn0,
            bgn1: self.len1 - self.end1,
            end1: self.len1 - self.bgn1,
            dist: self.dist,
            idt: self.idt,
            dist_c: self.dist_c,
            max_dist_c: self.max_dist_c,
            idt_c: self.idt_c,
            flag: self.flag | 0x80,
        }
    }
}

pub type OverlapMap = FxHashMap<u32, Vec<Overlap>>;

fn get_upath(g: &OvlpGraph, v: ReadNode, w: ReadNode) -> Vec<ReadNode> {
    //
    // find a simple path (no out branch) start at v->w edge in the graph g
    //

    let mut path = Vec::<ReadNode>::new();
    let mut visited_nodes = FxHashSet::<ReadNode>::default();
    path.push(v);
    path.push(w);

    visited_nodes.insert(v);
    visited_nodes.insert(w);

    let mut n = w;
    loop {
        if g.neighbors_directed(n, Outgoing).count() == 1
            && g.neighbors_directed(n, Incoming).count() == 1
        {
            n = *g
                .neighbors_directed(n, Outgoing)
                .into_iter()
                .collect::<Vec<ReadNode>>()
                .get(0)
                .unwrap();
            if !visited_nodes.contains(&n) {
                path.push(n);
                visited_nodes.insert(n);
            } else {
                path.push(n);
                break;
            }
        } else {
            break;
        }
    }
    path
}

pub fn get_utg_paths(g: &OvlpGraph) -> Vec<(u32, Vec<ReadNode>)> {

    //
    // get all unitig paths
    //

    let mut start_nodes = FxHashSet::<(u32, u8)>::default();
    for v in g.nodes() {
        if g.neighbors_directed(v, Incoming).count() != 1
            || g.neighbors_directed(v, Outgoing).count() != 1
        {
            start_nodes.insert(v);
        }
    }

    let mut uid = 0_u32;
    let mut paths = Vec::<(u32, Vec<(u32, u8)>)>::new();

    for v in start_nodes {
        for w in g.neighbors_directed(v, Outgoing) {
            let path = get_upath(&g, v, w);
            //let e = path[path.len() - 1];
            paths.push((uid, path));
            uid += 1;
        }
    }
    paths
}

pub fn transitive_reduction(g: &mut U32AsmGraph) -> () {    
    let mut tr_edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
    for v in g.nodes().into_iter() {
        let mut edges = Vec::<(u32, (u32, u8), (u32, u8))>::with_capacity(32);
        for w in g.neighbors_directed(v, Outgoing) {
            let ovlp_length = *g.edge_weight(v, w).unwrap();
            edges.push((ovlp_length, v, w));
        }
        edges.sort();
        for (_l, v, w) in edges.iter() {
            for (_l, _vv, x) in edges.iter() {
                if x == w {
                    continue;
                }
                // w in out_neighbor(x), v->x->w exist, v->w can be eliminate
                let mut found = false;
                for x_out in g.neighbors_directed(*x, Outgoing) {
                    if x_out == *w {
                        found = true;
                        break;
                    }
                }
                if found {
                    tr_edges.insert((*v, *w));
                }
            }
        }
    }

    for (v, w, _) in g.all_edges() {
        if g.neighbors_directed(v, Incoming).count() == 0
            && g.neighbors_directed(w, Incoming).count() >= 2
            && g.neighbors_directed(w, Outgoing).count() >= 1
        {
            tr_edges.insert((v, w));
        }
        if g.neighbors_directed(w, Outgoing).count() == 0
            && g.neighbors_directed(v, Outgoing).count() >= 2
            && g.neighbors_directed(v, Incoming).count() >= 1
        {
            tr_edges.insert((v, w));
        }
    }

    for (v, w) in tr_edges.into_iter() {
        //println!("D {:?} {:?}", v, w);
        g.remove_edge(v, w);
    }
}

pub fn remove_simple_spur(g: &mut U32AsmGraph, max: u32) {
    let mut edges = FxHashSet::<OvlpEdge>::default();
    for (v, w, s) in g.all_edges() {
        if *s > max {
            continue;
        }
        if g.neighbors_directed(v, Incoming).count() == 0
            && g.neighbors_directed(w, Incoming).count() >= 2
            && g.neighbors_directed(w, Outgoing).count() >= 1
        {
            edges.insert((v, w));
        }
        if g.neighbors_directed(w, Outgoing).count() == 0
            && g.neighbors_directed(v, Outgoing).count() >= 2
            && g.neighbors_directed(v, Incoming).count() >= 1
        {
            edges.insert((v, w));
        }
    }
    for (v, w) in edges.into_iter() {
        g.remove_edge(v, w);
    }
}

pub fn remove_single_bridge(g: &mut U32AsmGraph, max: u32) -> () {
    //
    // implement a huristic rule removing spurious connections (repeat induced bridge)
    //

    let mut to_remove = FxHashSet::<(u32, OvlpEdge)>::default();

    for (v, w, weight) in g.all_edges() {
        if g.neighbors_directed(v, Outgoing).count() >= 2
            && g.neighbors_directed(v, Incoming).count() >= 1
            && g.neighbors_directed(w, Outgoing).count() >= 1
            && g.neighbors_directed(w, Incoming).count() >= 2
            && *weight < max
        {
            to_remove.insert((*weight, (v, w)));
        }
    }

    let mut to_remove_vec = Vec::<(u32, OvlpEdge)>::new();
    for e in to_remove {
        to_remove_vec.push(e);
    }

    to_remove_vec.sort(); // remove edges with smaller weight (~ overlaps length) first

    for (_c, (v, w)) in to_remove_vec {
        if g.neighbors_directed(v, Outgoing).count() > 1
            && g.neighbors_directed(w, Incoming).count() > 1
            && g.neighbors_directed((w.0, 1 - w.1), Outgoing).count() > 1
            && g.neighbors_directed((v.0, 1 - v.1), Incoming).count() > 1
        {
            g.remove_edge(v, w);
            g.remove_edge((w.0, 1 - w.1), (v.0, 1 - v.1));
            //println!("RMSB0 {}:{} {}:{}", v.0, v.1, w.0, w.1);
            //println!("RMSB1 {}:{} {}:{}", w.0, 1-w.1, v.0, 1-v.1);
        }
    }
}

pub fn dump_utg_paths(
    paths: &Vec<(u32, Vec<ReadNode>)>,
    utg_g: &UtgGraph,
    filename: &String,
) -> Result<(), io::Error> {
    let mut utg_file = BufWriter::new(File::create(filename).unwrap());

    for (uid, p) in paths.iter() {
        let v = p[0];
        let w = p[p.len() - 1];
        let mut tag = 0_u32;
        if utg_g.contains_edge(v, w) {
            tag = 1;
        }
        writeln!(
            utg_file,
            "UTG {} {}:{} {}:{} {} {}",
            uid,
            v.0,
            v.1,
            w.0,
            w.1,
            p.len(),
            tag
        )?;
        for v in p {
            writeln!(utg_file, "N {} {}:{}", uid, v.0, v.1)?;
        }
    }
    Ok(())
}

pub fn dump_utg_gfa(
    paths: &Vec<(u32, Vec<ReadNode>)>,
    utg_g: &UtgGraph,
    rpair2overlap: &FxHashMap<ReadPair, Overlap>,
    _read_to_ctg: &FxHashMap<ReadNode, Vec<(u32, u32)>>,
    filename: &String,
) -> Result<(), io::Error> {
    let mut gfa_file = BufWriter::new(File::create(filename).unwrap());
    writeln!(gfa_file, "H\tVN:Z:1.0")?;
    let mut end_nodes = FxHashMap::<(u32, u8), Vec<u32>>::default();
    let mut bgn_nodes = FxHashMap::<(u32, u8), Vec<u32>>::default();
    let mut read_len = FxHashMap::<u32, u32>::default();
    for (uid, p) in paths.iter() {
        if p.len() < 2 {
            continue;
        }
        let nb = p[0];
        let ne = p[p.len() - 1];
        let mut tag = 0_u32;
        if utg_g.contains_edge(nb, ne) {
            tag = 1;
        }
        bgn_nodes.entry(nb).or_insert_with(|| vec![]).push(*uid);
        end_nodes.entry(ne).or_insert_with(|| vec![]).push(*uid);
        let mut node_start = Vec::<(ReadNode, u32)>::new();
        let mut v = p[0];
        let mut utg_len = 0_u32;
        node_start.push((v, 0));
        for w in p[1..p.len()].iter() {
            let ovlp = rpair2overlap.get(&ReadPair::new(v, *w)).unwrap();
            let len0 = ovlp.len0;
            let len1 = ovlp.len1;
            let ovlp_len = ovlp.end0 - ovlp.bgn0;
            if utg_len == 0 {
                utg_len += len0;
                read_len.insert(v.0, len0);
            }
            node_start.push((*w, utg_len));
            utg_len += len1 - ovlp_len;
            read_len.insert(w.0, len1);
            v = *w;
        }

        writeln!(
            gfa_file,
            "S\t{}\t*\tLN:i:{}\tRB:i:{}\tSB:i:{}\tRE:i:{}\tSE:i:{}\tNC:i:{}\tTG:i:{}",
            uid,
            utg_len,
            nb.0,
            nb.1,
            ne.0,
            ne.1,
            p.len(),
            tag
        )?;

        for (n, start) in node_start {
            writeln!(gfa_file, "N\t{}\t{}\t{}\t{}", uid, n.0, n.1, start)?;
        }
    }
    for (v, uids0) in end_nodes.iter() {
        if let Some(uids1) = bgn_nodes.get(v) {
            for u0 in uids0 {
                for u1 in uids1 {
                    let rlen = read_len.get(&v.0).unwrap();
                    writeln!(gfa_file, "L\t{}\t+\t{}\t+\t{}M", u0, u1, rlen)?;
                }
            }
        }

        if let Some(uids1) = end_nodes.get(&(v.0, 1 - v.1)) {
            for u0 in uids0 {
                for u1 in uids1 {
                    if *u0 == *u1 {
                        continue;
                    }
                    let rlen = read_len.get(&v.0).unwrap();
                    writeln!(gfa_file, "L\t{}\t+\t{}\t-\t{}M", u0, u1, rlen)?;
                }
            }
        }
    }
    Ok(())
}

pub fn bfs_extend(v: (u32, u8), g: &U32AsmGraph, limit: u32) -> (u32, Vec<(u32, u8)>) {
    //
    // output bfs search from node v in graph g bounded by the `limit`
    //

    let mut bfs = Bfs::new(g, v);
    let mut count = 0_u32;
    let mut nodes = Vec::<(u32, u8)>::with_capacity(32);
    while let Some(n) = bfs.next(g) {
        count += 1;
        if count >= limit {
            break;
        }
        nodes.push(n);
    }
    //assert!(nodes[0]==v);
    (count, nodes)
}

fn find_branching_nodes(
    g: &mut UtgGraph,
    max_path_length: u32,
    min_path_length: u32,
    max_edge_count: u32,
) -> Vec<(u32, u8)> {

    // output nodes that has branches

    let mut branching_nodes = FxHashSet::<(u32, u8)>::default();
    let candidates = g.nodes().collect::<Vec<ReadNode>>();
    for v in candidates {
        let out_count = g.neighbors_directed(v, Outgoing).count();
        if out_count < 2 {
            continue;
        }
        let mut ext_branch_count = 0_u32;
        let mut node_count = FxHashMap::<(u32, u8), u32>::default();
        for w in g.neighbors_directed(v, Outgoing) {
            let (count, p) = bfs_extend(w, g, max_path_length);
            let mut overlapped_path = false;
            let mut ovlp_count = 0;
            for ww in g.neighbors_directed(v, Outgoing) {
                ovlp_count += g.edge_weight(v, ww).unwrap();
            }
            for &vv in &p {
                for ww in g.neighbors_directed(vv, Outgoing) {
                    ovlp_count += g.edge_weight(vv, ww).unwrap();
                }
            }
            if ovlp_count <= 3 {
                // don't count if a branch is very short
                break;
            }
            if count >= min_path_length {
                let mut ec = 0_u32;
                for vv in p {
                    ec += 1;
                    if ec > max_edge_count {
                        break;
                    }
                    if vv.0 == v.0 {
                        // loop
                        break;
                    }
                    let c = node_count.entry(vv).or_insert(0);
                    *c += 1;
                    if *c > 1 {
                        overlapped_path = true;
                        break;
                    }
                }
            }
            if !overlapped_path {
                ext_branch_count += 1;
            }
        }
        if ext_branch_count >= 2 {
            branching_nodes.insert(v);
        }
    }
    branching_nodes.iter().map(|x| *x).collect()
}

fn find_path(g: &UtgGraph, s: ReadNode, t: ReadNode) -> Option<(u32, Vec<OvlpEdge>)> {

    // search a path from s to t

    assert!(s != t);
    let mut pre = FxHashMap::<ReadNode, Option<ReadNode>>::default();
    let mut min_edge_count = FxHashMap::<ReadNode, u32>::default();
    let mut out = Vec::<OvlpEdge>::new();

    pre.insert(s, None);
    min_edge_count.insert(s, 0);
    let limit = 256_u32;
    let mut count = 0_u32;
    let mut found = false;

    let mut bfs = Bfs::new(&g, s);
    let mut m: u32;
    while let Some(n) = bfs.next(&g) {
        m = u32::MAX;
        let mut min_in = None;
        for ww in g.neighbors_directed(n, Incoming) {
            if let Some(mm) = min_edge_count.get(&ww) {
                if *mm < m {
                    m = *mm;
                    min_in = Some(ww);
                }
            }
        }

        if min_in != None && n != s {
            pre.insert(n, min_in);
            min_edge_count.insert(n, m + g.edge_weight(min_in.unwrap(), n).unwrap());
        }

        if n == t {
            found = true;
            break;
        }
        count += 1;
        if count > limit {
            break;
        }
    }
    if found {
        let mut v = t;
        loop {
            let p = *pre.get(&v).unwrap();
            if p == None {
                break;
            }
            out.push((p.unwrap(), v));
            v = p.unwrap();
        }
        out.reverse();
        Some((*min_edge_count.get(&t).unwrap(), out))
    } else {
        None
    }
}

pub fn utg_reduction(paths: &Vec<(u32, Vec<ReadNode>)>, g0: &OvlpGraph) -> (UtgGraph, OvlpGraph) {

    // graph reduction in the unitig level

    let mut utg_g = UtgGraph::new();
    let mut g = OvlpGraph::new();
    for (_uid, p) in paths.iter() {
        let b = p[0];
        let e = p[p.len() - 1];
        utg_g.add_edge(b, e, p.len() as u32);
    }

    let mut branching_nodes = find_branching_nodes(&mut utg_g, 24, 1, 128);
    branching_nodes.sort();
    let mut bracnhing_node_sinks = FxHashSet::<ReadNode>::default();
    for v in branching_nodes.iter() {
        let w = (v.0, 1 - v.1);
        if utg_g.neighbors_directed(w, Outgoing).count() > 0 {
            // only add w if the w is not dead end
            bracnhing_node_sinks.insert(w);
            //println!("SINK {}:{}", v.0, 1 - v.1);
        }
    }
    //println!("BCOUNT {}", branching_nodes.len());
    let mut ripaths = Vec::<Vec<OvlpEdge>>::new(); //ripaths = repeat induced path
    for v in branching_nodes.iter() {
        let mut paths_to_remove_candidates = Vec::<(u32, Vec<OvlpEdge>)>::new();
        for w in utg_g.neighbors_directed(*v, Outgoing) {
            let vv = *v;
            let ww = w;
            if ww == vv {
                // ignore self edge caused by as loop
                continue;
            }
            log::debug!("BS {}:{} {}:{}", vv.0, vv.1, ww.0, ww.1);
            let c_nodes = bfs_extend(ww, &utg_g, 16);
            let mut sink_nodes = Vec::<ReadNode>::new();
            for p in &c_nodes.1 {
                if bracnhing_node_sinks.contains(p) {
                    sink_nodes.push(*p);
                }
            }
            if sink_nodes.len() > 0 {
                let mut s_ec = u32::MAX;
                let mut s_path = Vec::<OvlpEdge>::new();
                let ec0 = *utg_g.edge_weight(vv, ww).unwrap();
                s_path.push((vv, ww));
                for n in sink_nodes {
                    log::debug!("V {}:{} W {}:{}  SINK {}:{}", v.0, v.1, w.0, w.1, n.0, n.1);
                    if w == n {
                        s_ec = ec0;
                        break;
                    } else if let Some((mut ec, path)) = find_path(&utg_g, ww, n) {
                        ec += ec0;
                        if ec < 5 {
                            if ec < s_ec {
                                s_ec = ec;
                                s_path.extend(path);
                            }
                        }
                    }
                }
                if s_ec < 5 {
                    log::debug!("PL {} {}", s_ec, s_path.len());
                    paths_to_remove_candidates.push((s_ec, s_path));
                }
            }
        }
        //remove repeat candidates but ensure to keep one out
        paths_to_remove_candidates.sort();
        /*
        println!(
            "BN {}:{} {}",
            v.0,
            v.1,
            utg_g.neighbors_directed(*v, Outgoing).count()
        );
        */
        for (_ec, path) in paths_to_remove_candidates.iter() {
            log::debug!("RMPC {} {:?}", _ec, path,);
        }
        if paths_to_remove_candidates.len() < utg_g.neighbors_directed(*v, Outgoing).count() {
            for (_ec, path) in paths_to_remove_candidates.iter() {
                log::debug!("RMP {} {:?}", _ec, path);
                ripaths.push(path.clone());
            }
        } else {
            let l = paths_to_remove_candidates.len();
            if l >= 1 {
                for i in 0..l - 1 {
                    let (_ec, path) = paths_to_remove_candidates.get(i).unwrap();
                    log::debug!("RMP {} {:?}", _ec, path);
                    ripaths.push(path.clone());
                }
            }
        }
    }
    for p in ripaths {
        for (vv, ww) in p {
            utg_g.remove_edge(vv, ww);
            utg_g.remove_edge((ww.0, 1 - ww.1), (vv.0, 1 - vv.1));
            //log::debug!("RME {}:{} {}:{}", vv.0, vv.1, ww.0, ww.1);
        }
    }

    //patch some useful edge backs
    for (_uid, p) in paths.iter() {
        let b = p[0];
        let e = p[p.len() - 1];
        if utg_g.neighbors_directed(b, Outgoing).count() == 0
            && utg_g.neighbors_directed(b, Incoming).count() > 0
            && utg_g.neighbors_directed(e, Outgoing).count() > 0
            && utg_g.neighbors_directed(e, Incoming).count() == 0
        {
            utg_g.add_edge(b, e, p.len() as u32);
            utg_g.add_edge((e.0, 1 - e.1), (b.0, 1 - b.1), p.len() as u32);
            //log::debug!("ADD {}:{} {}:{}", b.0, b.1, e.0, e.1);
        }
    }

    //remove_simple_spur(&mut utg_g, 5);
    //remove_single_bridge(&mut utg_g, 5);

    for (_uid, p) in paths.iter() {
        let v = p[0];
        let w = p[p.len() - 1];
        if utg_g.contains_edge(v, w) {
            let mut vv = v;
            for i in 1..p.len() {
                let ww = p[i];
                let weight = *g0.edge_weight(vv, ww).unwrap();
                g.add_edge(vv, ww, weight);
                vv = ww;
            }
        }
    }
    (utg_g, g)
}
