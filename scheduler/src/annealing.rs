use rand;
use rand::distributions::{IndependentSample, Range};
use std;
use std::io::Write;
use std::thread;

pub trait Annealer {
    type State;
    type Move;

    fn init_state(&self) -> Self::State;
    fn start_temp(&self, init_score: f64) -> f64;

    fn is_done(&self, _score: f64) -> bool {
        false
    }

    fn eval(&self, state: &Self::State) -> f64;
    fn neighbour(&self, state: &Self::State) -> Self::Move;

    fn apply(&self, state: &mut Self::State, mov: &Self::Move);
    fn unapply(&self, state: &mut Self::State, mov: &Self::Move);

    fn apply_and_eval(&self, state: &mut Self::State, mov: &Self::Move, _prev_score: f64) -> f64 {
        self.apply(state, mov);
        self.eval(state)
    }
}

#[derive(Clone)]
pub struct AnnealingOptions {
    pub steps: usize,
    pub limit_temp: f64,
    pub restart: usize,
    pub threads: usize,
    pub silent: bool,
}

pub fn annealing<A: 'static + Annealer + Clone + Send + Sync>(
    annealer: &A,
    opt: &AnnealingOptions,
) -> (f64, <A as Annealer>::State)
where
    <A as Annealer>::State: Clone + Send + Sync,
{
    assert!(opt.threads > 0);

    if opt.threads == 1 {
        do_annealing(None, annealer, opt)
    } else {
        let mut ths = vec![];

        for i in 0..opt.threads {
            let a = annealer.clone();
            let o = opt.clone();
            ths.push(thread::spawn(move || do_annealing(Some(i), &a, &o)));
        }

        ths.into_iter()
            .map(|th| th.join().unwrap())
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
            .unwrap()
    }
}

fn do_annealing<A: Annealer>(
    thread_id: Option<usize>,
    annealer: &A,
    opt: &AnnealingOptions,
) -> (f64, <A as Annealer>::State)
where
    A::State: Clone,
{
    let mut rng = rand::thread_rng();
    let drng = Range::new(0.0, 1.0);

    let mut state = annealer.init_state();
    let mut cur_score = annealer.eval(&state);
    let mut best_score = cur_score;
    let mut best_ans = state.clone();

    macro_rules! progress {
        ($($arg:expr),*) => {
            if !opt.silent {
                let err = std::io::stderr();
                let mut err = err.lock();
                if let Some(tid) = thread_id {
                    write!(err, "[{:02}] ", tid).unwrap();
                }
                writeln!(err, $($arg),*).unwrap();
            }
        };
    }

    progress!("Initial score: {}", cur_score);

    let mut restart_cnt = 0;

    let t_max = annealer.start_temp(cur_score);
    let t_min = opt.limit_temp;

    let step = opt.steps as f64;
    let decay = ((t_min / t_max).ln() / step).exp();

    progress!("Temperature decay: {}", decay);

    let mut temp = t_max;
    let mut cnt = 0;
    loop {
        if temp < t_min {
            restart_cnt += 1;
            if restart_cnt >= opt.restart {
                break;
            }
            progress!("Restarting... {}/{}", restart_cnt, opt.restart);
            temp = 100.0;
        }

        cnt += 1;
        if cnt % 100000 == 0 {
            progress!("Current state: score = {:.3}, temp = {:.9}", cur_score, temp);
        }

        let mov = annealer.neighbour(&state);
        let new_score = annealer.apply_and_eval(&mut state, &mov, cur_score);

        if new_score <= cur_score
            || drng.ind_sample(&mut rng) <= ((cur_score - new_score) as f64 / temp).exp()
        {
            cur_score = new_score;
            if cur_score < best_score {
                if best_score - cur_score > 1e-6 {
                    progress!("Best: score = {:.3}, temp = {:.9}", cur_score, temp);
                }
                best_score = cur_score;
                best_ans = state.clone();
            }
            if annealer.is_done(cur_score) {
                break;
            }
        } else {
            annealer.unapply(&mut state, &mov);
        }

        temp *= decay;
    }
    (best_score, best_ans)
}
