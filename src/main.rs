
//uses https://gregstrohman.com/wp-content/uploads/2019/10/jic01-Interval-tuning.pdf as guideline for interval-based JI tuning.

//TODO: after getting REALLY SIMPLE midi export going for test purposes,
//make sure this tool interfaces well with the leading JUCE alternative in the ecosystem
//(that is, nih-plug). In particular, it should remain as agnostic as possible WRT data types, etc.

use nalgebra::{DMatrix, DVector};
use midly::Smf;

static ET_TO_JUST: [f64; 12] = [
    0.0, // 0, unis.
    1.1173, // 1
    2.0391, // 2, != 12 - ET_TO_JUST[10]
    3.1564, // 3
    3.8631, // 4
    4.9804, // 5
    5.8251, // 6, != 12 - ET_TO_JUST[6]
    //Tritone is subject to empirical revision: might just want to be a perfect 6.0!
    7.0196, // 7
    8.1369, // 8
    8.8436, // 9
    9.6883, // 10, != 12 - ET_TO_JUST[2]
    10.8827 // 11
]; //Maps integer intervals mod 12 to just intervals

static IMPORTANCE_WEIGHTS: [f64; 12] = [
    1.0, //unis
    0.25, //m2
    0.5, //M2
    1.0, //m3
    1.0, //M3
    1.0, //p4
    0.125, //tritone
    1.0, //p5
    1.0, //m6
    1.0, //M6
    0.25, //m7
    0.25, //M7
]; //Made-up weights, based on how musically important I think it is that these intervals sound in-tune.
//Another thought: weights could possibly be functions of octave distance?

fn weight_vec(equal_intervals: &Vec<i8>) -> Vec<f64> {
    let m = equal_intervals.len();
    let mut out = vec![];
    for i in 0..m {
        out.push(IMPORTANCE_WEIGHTS[(equal_intervals[i] % 12) as usize]);
    }
    out   
}

fn compute_tuning_vector_f64(equal_notes: &Vec<i8>) -> Result<Vec<f64>, String> {

    //Presumes the notes vector is sorted.

    // Least squares: choose tuning vector that minimizes
    // ||W * A * (tuning + equal) - W * just||^2,
    // Where A takes a vector of pitches to a vector holding
    // the interval from pitch 1 to 2, then 1 to 3, ... then 2 to 3, ... then n-1 to n.
    // that is, A equals:
    // [-1 1 0 0 ... 0]
    // [-1 0 1 0 ... 0]
    // [      ...     ]
    // [0 -1 1 0 ... 0]
    // [0 -1 0 1 ... 0]
    // [      ...     ]
    // [0 ... 0 0 -1 1]

    //W, conversely, is the diagonal matrix holding the square roots of the weights for each interval.

    //So, (tuning + equal) = (WA)^(+) * W * just solves the system.

    //Since the nullspace of A^T * A is span([1; 1; ... 1]), the best answer
    //will be that given by the above, plus some scalar multiple of [1; 1; ... 1].

    //Instead of the theoretically perfect solution (identifying the root somehow, and tuning the root to its 12tet pitch),
    //we instead choose the linear shift that involves the smallest mean square deviation from equal temperment across all notes.
    //so, choose a to minimize ||equal + (tuning + a * 1_vec) - equal||^2,
    //i.e., minimize ||1_vec * a + tuning||^2
    //This turns out to be equivalent to mean-centering tuning.

    let n = equal_notes.len();
    if n == 0 {return Ok(vec![]);}
    if n == 1 {return Ok(vec![0.0]);}
    let m = (n-1) * n / 2;
    
    let (
            just_intervals_weighted, //W * j
            WA //W * A
        ) = {
            let mut j_weighted= vec![];
            let mut WA = DMatrix::zeros(m, n);
            let mut start = 0;
            let mut end = 1;
            let mut i = 0;
            loop {
                let interval = (equal_notes[end] - equal_notes[start]) as usize;
                let current_weight = f64::sqrt(IMPORTANCE_WEIGHTS[interval % 12]);
                j_weighted.push(current_weight * ((interval - (interval % 12)) as f64 + ET_TO_JUST[interval % 12]));
                WA[(i, start)] = -current_weight;
                WA[(i, end)] = current_weight;
                end += 1;
                i += 1;
                if end == n {
                    start += 1;
                    end = start + 1;
                }
                if end == n {
                    break;
                }
            }
            (
                DVector::from_vec(j_weighted),
                WA
            )
        }; 

    //WA will be singular: SVD is the best bet for the pseudoinverse.
    let epsilon = 0.00001;
    let wa_svd = WA.svd(true, true);
    let svd_failed_msg = &String::from("Svd failed.");
    let u = wa_svd.u.ok_or(svd_failed_msg)?;
    let v_t = wa_svd.v_t.ok_or(svd_failed_msg)?;
    let mut s = DMatrix::from_diagonal(&wa_svd.singular_values);
    
    for i in 0..s.ncols() {
        if s[(i,i)] < epsilon {
            s[(i, i)] = 0.0;
        } else {
            s[(i, i)]  = 1.0 / s[(i, i)];
        }
    }

    let mut tuning =  v_t.transpose() * (s * (u.transpose() * just_intervals_weighted));

    let mut tuning_sum = 0.0;
    for i in 0..n {
        tuning[(i, 0)] -= equal_notes[i] as f64;
        tuning_sum += tuning[(i,0)];
    }

    let tuning_mean = tuning_sum / (n as f64);

    for i in 0..n {
        tuning[(i,0)] -= tuning_mean;
    }

    let mut tuning_out = vec![];

    for i in 0..n {
        tuning_out.push(tuning[(i,0)]);
    }

    Ok(tuning_out)

}

//Note: Midly has a function to convert floats to pitchbends (which are 14-bit ints???)

#[cfg(test)]
mod tests {
    //NOTE: rust runs tests in parallel! Use print statements at your peril!
    use super::*;
    const TEST_LO: i8 = 36; //C2
    const TEST_HI: i8 = 100; //E7
    const EPSILON: f64 = 0.0001;
    fn jury_rigged_random(a: i8) -> i8 {
        if rand::random() && a < 4 {
            jury_rigged_random(a + 1)
        } else {
            a
        }
    }
    //#[test]
    fn unison_tests() {
        //ensure tuning_vector for single notes is EMPTY (not just zero)
        // println!("Unison:");
        // for i in TEST_LO..TEST_HI {
        //     let equal_notes = vec![i];
        //     let tuning_vector = compute_tuning_vector_f64(&equal_notes, 0).unwrap();
        //     for j in 0..tuning_vector.len() {
        //         println!("tuning_vector[{}] = {}", j, tuning_vector[j]);
        //         assert_eq!(0.0, tuning_vector[j]);
        //     }
        // }
        //ensure tuning vector for octaves is zero
        println!("Octaves:");
        for i in TEST_LO..(TEST_HI / 2) {
            let r1 = jury_rigged_random(0);
            let r2 = jury_rigged_random(1);
            let equal_notes = 
                {let mut a = vec![i, (i + 12 * r1), (i + 12 * r2)];
                a.sort_unstable(); a};
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("equal_notes = [{}, {}, {}]", equal_notes[0], equal_notes[1], equal_notes[2]);
            println!("tuning_vector = [{}, {}, {}]", tuning_vector[0], tuning_vector[1], tuning_vector[2]);    
            for j in 0..tuning_vector.len() {
                assert!(tuning_vector[j]- 0.0 < EPSILON);
            }
        }
    }

    //#[test]
    fn interval_tests(){
        //ensure tuning vector for dyads is as given by the table.
        println!("Intervals:");
        for root in TEST_LO..TEST_HI- 12 {
            for interval in 0i8..27 {
                let top = root + interval;
                let equal_notes = vec![root, top];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                for j in 0..tuning_vector.len() {
                    assert!((tuning_vector[j] - ET_TO_JUST[(interval % 12) as usize]) < EPSILON)
                }
                println!("equal_notes = [{}, {}]", equal_notes[0], equal_notes[1]);
                println!("tuning_vector = [{}, {}]", tuning_vector[0], tuning_vector[1]);
            }
        }
    }
    #[test]
    fn triad_tests() {
        //requires auditory checks as well
        //major triad:
        //in root:
       println!("Triads -- root posision:");
       for root in 10..11 {
            let equal_notes = vec![root, root + 4, root + 7];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("equal_notes = [{}, {}, {}]", equal_notes[0], equal_notes[1], equal_notes[2]);
            println!("tuning_vector = [{}, {}, {}]", tuning_vector[0], tuning_vector[1], tuning_vector[2]);
            //play audio
       }
       //in 1st:
       println!("Triads -- first inversion:");
       for root in 10..11 {
            let equal_notes = vec![root + 4, root + 7, root + 12];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("equal_notes = [{}, {}, {}]", equal_notes[0], equal_notes[1], equal_notes[2]);
            println!("tuning_vector = [{}, {}, {}]", tuning_vector[0], tuning_vector[1], tuning_vector[2]);
            //play audio
       }
       //in 2nd:
       println!("Triads -- second inversion:");
       for root in 10..11 {
            let equal_notes = vec![root + 7, root + 12, root + 16];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("equal_notes = [{}, {}, {}]", equal_notes[0], equal_notes[1], equal_notes[2]);
            println!("tuning_vector = [{}, {}, {}]", tuning_vector[0], tuning_vector[1], tuning_vector[2]);
            //play audio
       }
       
    }

    //#[test]
    fn stack_tests() {
        //p4 x5
        let max = 3;
        let mut equal_notes = vec![0i8];
        for i in 0..max {
            equal_notes.push(equal_notes[i] + 5);
        }
        let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
        println!("Stacked 4ths: notes: {:?}", equal_notes);
        println!("tuning: {:?}", tuning_vector);

        //p5 x5
        let mut equal_notes = vec![0i8];
        for i in 0..max {
            equal_notes.push(equal_notes[i] + 7);
        }
        let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
        println!("stacked 5ths: notes: {:?}", equal_notes);
        println!("tuning: {:?}", tuning_vector);
    }

    //non-rigorous-- sanity check to see how the algorithm affects more complicated chords
    #[test]
    fn extension_tests() {
        //dom7
        let equal_notes = vec![0i8, 4, 7, 10];
        let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
        println!("dominant 7: notes: {:?}", equal_notes);
        println!("tuning: {:?}", tuning_vector);
        //lowers the fifth by a ton: not ideal!

        //m9
        let equal_notes = vec![0i8, 7, 10, 14, 15];
        let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
        println!("minor 9: notes: {:?}", equal_notes);
        println!("tuning: {:?}", tuning_vector);
        //seems normal

        //M7#11
        let equal_notes = vec![0i8, 7, 11, 16, 18];
        let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
        println!("maj7#11: notes: {:?}", equal_notes);
        println!("tuning: {:?}", tuning_vector);
        //seems normal
    }

    fn asymmetry_tests() {
       //test on things like dominant 7th chords with 2 roots represented, to ensure it doesn't drag up the 7th too much.
       
    }
    
    
}

fn main() {

    use std::env;
    use std::fs;

    let args:Vec<String> = env::args().collect();

    let (midi_bytes, midi_names) = {
        let mut out = vec![];
        let mut names = vec![];
        for arg in 1..args.len() {
            if let Ok(next) = fs::read(&args[arg]) {
                out.push(next);
                names.push(&args[arg]);
            } else {
                panic!("Failed to read file \'{}\'.", &args[arg]);
            }
        }
        (out, names)
    };

    let Smf_files = {
        let mut out = vec![];
        for i in 0..midi_bytes.len() {
            if let Ok(next) = Smf::parse(&midi_bytes[i]) {
                out.push(next);
            } else  {
                panic!("Failed to parse file \'{}\' as midi.", midi_names[i]);
            }
        }
        out
    };

    //midly is the rust midi parser/writer you want to use!!

}