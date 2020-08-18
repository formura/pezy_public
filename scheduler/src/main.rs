extern crate rand;
extern crate regex;
#[macro_use]
extern crate structopt;

use rand::thread_rng;
use rand::distributions::{IndependentSample, Range};
use regex::Regex;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::path::PathBuf;
use structopt::StructOpt;

mod annealing;
use annealing::*;

type Result<T> = std::result::Result<T, Box<std::error::Error>>;

/*

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "lowercase", tag = "type")]
enum Inst {
    Alu {
        loc: usize,
        opc: String,
        args: Vec<usize>,
    },
    I32 {
        loc: usize,
        value: i32,
    },
    F64 {
        loc: usize,
        value: f64,
    },
    Var {
        loc: usize,
        ident: String,
    },
    Load {
        loc: usize,
        ident: String,
        indices: Vec<usize>,
    },
    Store {
        loc: usize,
        ident: String,
        indices: Vec<usize>,
    }
}

impl Inst {
    fn loc(&self) -> usize {
        use Inst::*;
        match self {
            &Alu { loc, .. } => loc,
            &I32 { loc, .. } => loc,
            &F64 { loc, .. } => loc,
            &Var { loc, .. } => loc,
            &Load { loc, .. } => loc,
            &Store { loc, .. } => loc,
        }
    }

    fn refs(&self) -> &[usize] {
        use Inst::*;
        match self {
            &Alu { ref args, .. } => &args,
            &Load { ref indices, .. } => &indices,
            &Store { ref indices, .. } => &indices,
            _ => &[],
        }
    }
}

fn liveness_anal(insts: &[Inst], mapping: &[usize], penalty: bool) -> f64 {
    let mut ret = 0;
    // let mut ss = BTreeSet::new();
    let mut ss = vec![false; insts.len() + 256];
    let mut hist = vec![0; insts.len() + 256];
    let mut cur = 0;
    let mut cost = 0.0;

    for ii in (0..insts.len()).rev() {
        let i = mapping[ii];

        for r in insts[i].refs().iter() {
            if !ss[*r] {
                ss[*r] = true;
                cur += 1;
            }
        }

        let loc = insts[i].loc();
        if ss[loc] {
            ss[loc] = false;
            cur -= 1;
        }

        ret = ret.max(cur);
    }

    if penalty {
        for ii in 0..insts.len() {
            let i = mapping[ii];

            for r in insts[i].refs().iter() {
                cost += (ii - hist[*r]).pow(2) as f64;
            }

            let loc = insts[i].loc();
            hist[loc] = ii;
        }
    }

    ret as f64 + cost / hist.len().pow(2) as f64
}

fn can_swap(a: &Inst, b: &Inst) -> bool {
    !b.refs().contains(&a.loc())
}

struct MultiSet<T>(BTreeSet<T, u32>);

impl<T: Ord> MultiSet<T> {

}

struct AnalState {
    mapping: Vec<usize>,
    lives: Vec<Vec<usize>>,
    live_cnts: BTreeMap<i32, i32>;
}

impl AnalState {
    fn new(insts: &[Inst], mapping: &[usize]) -> AnalState {
        let mut ss = vec![0; insts.len() + 256];
        let mut cur = 0;

        let mut lives = vec![];
        let mut live_cnts = vec![];

        for ii in (0..insts.len()).rev() {
            let i = mapping[ii];

            for r in insts[i].refs().iter() {
                ss[*r] += 1;
                if ss[*r] == 1{
                    cur += 1;
                }
            }

            let loc = insts[i].loc();
            if ss[loc] > 0 {
                ss[loc] = 0;
                cur -= 1;
            }

            lives.push(ss.clone());
            live_cnts.push(cur);
        }

        lives.rev();
        live_cnts.rev();

        AnalState {
            mapping: mapping.to_owned(),
            lives,
            live_cnts,
        }
    }

    fn swap(&mut self, insts: &[Inst], a: usize, b: usize) -> i32 {
        assert!(a + 1 == b);

        let ma = self.mapping[a];
        let mb = self.mapping[b];

        let mut diff = 0;

        diff += self.unapply(a, &insts[ma]);
        // diff += self.unapply(a, &insts[mb]);
        // diff += self.unapply(b, &insts[mb]);

        // diff += self.apply(a, &insts[ma]);
        // diff += self.apply(a, &insts[mb]);
        diff += self.apply(b, &insts[ma]);

        self.mapping.swap(a, b);
    }

    fn unapply(&mut self, i: usize, inst: &Inst) -> i32 {
        for r in inst.refs().iter() {
            ss[*r] += 1;
            if ss[*r] == 1{
                cur += 1;
            }
        }

        let loc = insts[i].loc();
        if ss[loc] > 0 {
            ss[loc] = 0;
            cur -= 1;
        }

    }
}
*/

#[derive(Debug, Clone)]
enum Line {
    Inst {
        opc: String,
        args: Vec<String>,
        srcs: Vec<String>,
        dest: String,
    },
    Label {
        label: String,
    }
}

fn to_res(s: &str) -> &str {
    if s.chars().nth(0).unwrap() == 'x' || s.chars().nth(0).unwrap() == 'd' {
        return s;
    }

    if s.chars().nth(0).unwrap().is_digit(10) {
        return "imm"
    }

    match s {
        "sp" => "x2",
        "zr" => "x0",

        "r0" | "r1" => "x0",
        "r2" | "r3" => "x1",
        "r4" | "r5" => "x2",
        "r6" | "r7" => "x3",
        "r8" | "r9" => "x4",
        "r10" | "r11" => "x5",
        "r12" | "r13" => "x6",
        "r14" | "r15" => "x7",
        "r16" | "r17" => "x8",
        "r18" | "r19" => "x9",
        "r20" | "r21" => "x10",
        "r22" | "r23" => "x11",
        "r24" | "r25" => "x12",
        "r26" | "r27" => "x13",
        "r28" | "r29" => "x14",
        "r30" | "r31" => "x15",

        _ => unreachable!("{}", s),
    }
}

impl Line {
    fn mk_inst(opc: String, args: Vec<String>) -> Line {
        let t = Line::Inst {
            opc: opc.clone(),
            args: args.clone(),
            srcs: vec![],
            dest: String::default(),
        };

        let srcs: Vec<String> = t.srcs().into_iter().map(|s| s.to_owned()).collect();
        let dest = if opc == "c.bf" {
            "NA".to_owned()
        } else {
            t.dst().to_owned()
        };

        Line::Inst {
            opc,
            args,
            srcs,
            dest,
        }
    }

    fn srcs(&self) -> Vec<&str> {
        if let &Line::Inst { ref opc, ref args, .. } = self {
            if opc == "d.ldd" || opc == "d.eldd" {
                let mut t: Vec<&str> = args[1..].iter().map(|s| to_res(s)).collect();
                t.push("mem");
                t
            } else if args.len() >= 1 {
                args[1..].iter().map(|s| to_res(s)).collect()
            } else {
                vec![]
            }
        } else {
            unreachable!()
        }
    }

    fn dst(&self) -> &str {
        if let &Line::Inst { ref opc, ref args, .. } = self {
            if opc == "d.sd" || opc == "d.esd" {
                "mem"
            } else if args.len() >= 1 {
                to_res(&args[0])
            } else {
                "***"
            }
        } else {
            unreachable!()
        }
    }

    fn latency(&self) -> i32 {
        if let &Line::Inst { ref opc, .. } = self {
            match opc.as_ref() {
                "d.mad" | "d.nmad" | "d.msub" => 2,
                "d.ldd" | "d.sd" => 4,
                _ => 1,
            }
        } else {
            unreachable!()
        }
    }

    fn to_string(&self) -> String {
        match self {
            &Line::Inst { ref opc, ref args, .. } => {
                format!("\t{}\t{}", opc, args.join(" "))
            }
            &Line::Label { ref label } => {
                format!("{}", label)
            }
        }
    }
}

fn is_mad(opc: &str) -> bool {
    opc == "d.mad" || opc == "d.nmad" || opc == "d.nmad"
}

fn is_mem(opc: &str) -> bool {
    opc == "d.ldd" || opc == "d.eldd" || opc == "d.sd" || opc == "d.esd"
}

fn can_swap(l1: &Line, l2: &Line) -> bool {
    use Line::*;

    match (l1, l2) {
        (&Inst {opc: ref opc1, srcs: ref s1, dest: ref d1, .. },
         &Inst {opc: ref opc2, srcs: ref s2, dest: ref d2, .. }) => {
            if opc1 == "c.bf" || opc2 == "c.bf" {
                return false;
            }

            !s2.contains(&d1) && !s1.contains(&d2)
        }
        _ => false
    }
}

fn can_dual(l1: &Line, l2: &Line) -> bool {
    use Line::*;

    if !can_swap(l1, l2) {
        return false;
    }

    // assert!(can_swap(l2, l1), "{:?} <=> {:?}", l2, l1);

    match (l1, l2) {
        (&Inst { opc: ref opc1, args: ref args1, .. }, &Inst { opc: ref opc2, args: ref args2, .. }) => {
            if is_mad(opc1) {
                return is_mem(opc2)
            }
            if is_mad(opc2) {
                return is_mem(opc1)
            }

            if is_mem(opc1) && is_mem(opc2) {
                return false;
            }

            true
        }
        _ => unreachable!()
    }
}

fn calc_score(insts: &[Line], perm: &[usize]) -> i32 {
    let mut ret = 0;
    let mut skip = false;

    let mut dual_cnt = 0;

    for i in 0..insts.len() {
        if skip {
            skip = false;
            continue;
        }

        if i + 1 >= insts.len() {
            ret += 1;
            continue;
        }

        let op1 = &insts[perm[i]];
        let op2 = &insts[perm[i + 1]];

        if can_dual(op1, op2) {
            ret += 1;
            dual_cnt += 1;
            skip = true;
            // eprintln!("D: {:?}\n   {:?}", insts[perm[i]], insts[perm[i+1]]);
        } else {
            ret += 1;
            // eprintln!("S: {:?}", insts[perm[i]]);
        }
    }

    // eprintln!("Total instructions: {}", insts.len());
    // eprintln!("Dual issue cnt:     {}", dual_cnt);

    ret
}

fn parse_line(l: &str) -> Option<Line> {
    let l = if let Some(ix) = l.find("//") {
        l[0..ix].trim()
    } else {
        l
    };

    if l == "" {
        return None;
    }

    let mut it = l.split_whitespace();

    let fst = it.next().unwrap().to_owned();

    if fst.chars().rev().nth(0).unwrap() == ':' {
        return Some(Line::Label {
            label: fst
        });
    }

    let args = it.map(|s| s.to_owned()).collect();

    Some(Line::mk_inst(fst, args))
}

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(name = "FILE", parse(from_os_str))]
    file: PathBuf,
}

#[derive(Clone)]
struct Problem {
    insts: Vec<Line>,
}

impl Annealer for Problem {
    type State = Vec<usize>;
    type Move = (usize, usize);

    fn init_state(&self) -> Self::State {
        (0..self.insts.len()).collect()
    }

    fn start_temp(&self, init_score: f64) -> f64 {
        5.0
    }

    fn eval(&self, state: &Self::State) -> f64 {
        calc_score(&self.insts, &state) as f64
    }

    fn neighbour(&self, state: &Self::State) -> Self::Move {
        loop {
            let a = Range::new(0, state.len() - 1).ind_sample(&mut thread_rng());
            let b = a + 1;

            if can_swap(&self.insts[state[a]], &self.insts[state[b]]) {
                return (a, b)
            }
        }
    }

    fn apply(&self, state: &mut Self::State, mov: &Self::Move) {
        state.swap(mov.0, mov.1)
    }

    fn unapply(&self, state: &mut Self::State, mov: &Self::Move) {
        state.swap(mov.0, mov.1)
    }

    // fn apply_and_eval(&self, state: &mut Self::State, mov: &Self::Move, _prev_score: f64) -> f64 {
    //     self.apply(state, mov);
    //     self.eval(state)
    // }
}

fn app_main() -> Result<()> {
    let opt = Opt::from_args();

    let mut file = File::open(&opt.file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    let re = Regex::new("// step start\n((.|\\s)+)// step end")?;
    let caps = re.captures(&contents).unwrap();
    let cap = caps.get(0).unwrap();

    let asms = cap.as_str().to_owned();
    let lines: Vec<Line> = asms.lines().filter_map(|s| parse_line(s.trim())).collect();

    // println!("{:?}", lines);

    // println!("score: {}", calc_score(&lines, &(0..lines.len()).collect::<Vec<_>>()));

    println!("Start annealing...");

    let problem = Problem {
        insts: lines.clone(),
    };

    let (score, state) = annealing(&problem, &AnnealingOptions {
        steps: 100000,
        limit_temp: 1e-6,
        restart: 1,
        threads: 1,
        silent: false,
    });

    // println!("{:?}", state);
    // println!("Final score: {}", calc_score(&state));

    let mut ss = vec![];
    for &ix in state.iter() {
        ss.push(lines[ix].to_string());
    }

    let new_contents = re.replace(&contents,
        &("// step start optimized\n".to_string() + &ss.join("\n") + "\n// step end") as &str).to_owned();
    {
        let mut t = opt.file.clone();
        t.set_extension("opt.s");
        let mut file = File::create(&t)?;
        writeln!(file, "{}", new_contents);
    }

    Ok(())
}

fn main() {
    app_main().unwrap();
}
