
//uses https://gregstrohman.com/wp-content/uploads/2019/10/jic01-Interval-tuning.pdf as guideline for interval-based JI tuning.

use once_cell::sync::Lazy;
use nalgebra::{DMatrix, DVector};

static ET_TO_JUST: Lazy<[f64; 12]> = Lazy::new(|| {
    let mut arr = [0.0; 12];
    arr[0] = 0.0;
    arr[1] = 1.1173; 
    arr[2] = 2.0391; // != 12 - ET_TO_JUST[10]
    arr[3] = 3.1564;
    arr[4] = 3.8631;
    arr[5] = 4.9804;
    arr[6] = 5.8251; // != 12 - ET_TO_JUST[6]
    //Tritone is subject to empirical revision: might just want to be a perfect 6.0!
    arr[7] = 7.0196;
    arr[8] = 8.1369;
    arr[9] = 8.8436;
    arr[10] = 9.6883; // != 12 - ET_TO_JUST[2]
    arr[11] = 10.8827;
    arr
}); //Maps integer intervals mod 12 to just intervals

//Importance weights are never used without sqrting them first,
//but we leave them in the form of a comment for reference.

// static IMPORTANCE_WEIGHTS: [f64; 12] = [
//     1.0, //unis
//     0.25, //m2
//     0.5, //M2
//     1.0, //m3
//     1.0, //M3
//     1.0, //p4
//     0.125, //tritone
//     1.0, //p5
//     1.0, //m6
//     1.0, //M6
//     0.25, //m7
//     0.25, //M7
// ]; //Made-up weights, based on how musically important I think it is that these intervals sound in-tune.
// //Another thought: weights could possibly be functions of octave distance?

//An interesting thought: in a real soft synth, we could easily make these weights user-specified!

static IMPORTANCE_WEIGHTS_SQRTED: Lazy<[f64; 12]> = Lazy::new(|| {
    let mut arr = [0.0; 12];
    arr[0] = f64::sqrt(1.0); //unis
    arr[1] = f64::sqrt(0.25); //m2
    arr[2] = f64::sqrt(0.5); //M2
    arr[3] = f64::sqrt(1.0); //m3
    arr[4] = f64::sqrt(1.0); //M3
    arr[5] = f64::sqrt(1.0); //p4
    arr[6] = f64::sqrt(0.125); //tritone
    arr[7] = f64::sqrt(1.0); //p5
    arr[8] = f64::sqrt(1.0); //m6
    arr[9] = f64::sqrt(1.0); //M6
    arr[10] = f64::sqrt(0.5); //m7
    arr[11] = f64::sqrt(0.5); //M7
    arr
});

fn compute_tuning_vector_f64(equal_notes: &Vec<u8>) -> Result<Vec<f64>, String> {

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
                let current_weight = IMPORTANCE_WEIGHTS_SQRTED[interval % 12];
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

    //WA will be singular: SVD is a safe bet for the pseudoinverse.

    let wa_svd = WA.svd(true, true);
    let svd_failed_msg = &String::from("Svd failed.");
    let u = wa_svd.u.ok_or(svd_failed_msg)?;
    let v_t = wa_svd.v_t.ok_or(svd_failed_msg)?;
    let mut s = DMatrix::from_diagonal(&wa_svd.singular_values);
    let epsilon = 0.000001;
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
    use super::*;

    //MIDI convention: C4 is 60, etc.
    const TEST_LO: u8 = 21; //A0
    const TEST_HI: u8 = 100; //E7
    const EPSILON: f64 = 0.0001;
    fn within_epsilon_of(a: f64, b: f64) -> bool {
        return f64::abs(a - b) < EPSILON;
    }

    fn jury_rigged_random(a: u8) -> u8 {
        if rand::random() && a < 4 {
            jury_rigged_random(a + 1)
        } else {
            a
        }
    }
    mod numerical_tests {
        use super::*;
        #[test]
        fn unison_tests() {
            //ensure tuning vector for octaves and unison is the zero vector
            for i in TEST_LO..(TEST_HI / 2) {
                let r1 = jury_rigged_random(0);
                let r2 = jury_rigged_random(1);
                let equal_notes = 
                    {let mut a = vec![i, (i + 12 * r1), (i + 12 * r2)];
                    a.sort_unstable(); a};
                let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                for j in 0..tuning_vector.len() {
                    assert!(tuning_vector[j] < EPSILON);
                }
            }
        }

        #[test]
        fn interval_tests(){
            //ensure tuning vector for dyads is as given by the table,
            //after accounting for the vertical shift.
            //Also tests across multi-octave gaps.
            let mut failures = String::new();
            for root in TEST_LO..(TEST_HI - 12) {
                for interval in 0u8..27 {
                    let top = root + interval;
                    let equal_notes = vec![root, top];
                    let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                    let tuned_interval = (tuning_vector[1] + equal_notes[1] as f64) - (tuning_vector[0] + equal_notes[0] as f64);
                    let expected_interval = ET_TO_JUST[(interval % 12) as usize] + (12 * (interval / 12)) as f64;
                    if !within_epsilon_of(tuned_interval, expected_interval) {
                        failures += (format!("tuned result is {}, but expected {}. \n",
                        tuned_interval, expected_interval)).as_str();
                    }
                }
            }
            assert!(failures.is_empty(), "{failures}");
        }

        #[test]
        fn triad_tests() {
            //In real-world performance, major triads are typically mapped such that the M3
            //and the P5 are exactly as given by the tuning chart.
            //In this test, we check if our algorithm replicates this behavior.
            //That is, we check to see if our tuned vector is a vertical shift away from 
            //a major triad tuned in this way.
            //Root position:
            let mut failures = String::new();
            for root in 10..40 {
                let equal_notes = vec![root, root + 4, root + 7];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                let just_notes = {
                    let mut out = vec![];
                    for i in 0..3 {
                        out.push(equal_notes[i] as f64 + tuning_vector[i]);
                    }
                    out
                };
                if !within_epsilon_of(just_notes[1] - just_notes[0], ET_TO_JUST[4]) {
                    failures += format!("Root position triad failed to tune perfectly: third was {}, but expected {}", 
                    just_notes[1] - just_notes[0], ET_TO_JUST[4]).as_str();
                }
                if !within_epsilon_of(just_notes[2] - just_notes[0], ET_TO_JUST[7]) {
                    failures += format!("Root position triad failed to tune perfectly: fifth was {}, but expected {}", 
                    just_notes[1] - just_notes[0], ET_TO_JUST[4]).as_str();
                }
            }
            //1st inversion:
            for root in 10..40 {
                let equal_notes = vec![root + 4, root + 7, root + 12];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                let just_notes = {
                    let mut out = vec![];
                    for i in 0..3 {
                        out.push(equal_notes[i] as f64 + tuning_vector[i]);
                    }
                    out
                };
                if !within_epsilon_of(just_notes[2] - just_notes[0], ET_TO_JUST[8]) {
                    failures += format!("First inversion triad failed to tune perfectly: m6 was {}, but expected {}", 
                    just_notes[2] - just_notes[0], ET_TO_JUST[8]).as_str();
                }
                if !within_epsilon_of(just_notes[2] - just_notes[1], ET_TO_JUST[5]) {
                    failures += format!("First inversion triad failed to tune perfectly: p4 was {}, but expected {}", 
                    just_notes[1] - just_notes[0], ET_TO_JUST[4]).as_str();
                }
            }
            //2nd inversion:
            for root in 10..40 {
                let equal_notes = vec![root + 7, root + 12, root + 16];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
                let just_notes = {
                    let mut out = vec![];
                    for i in 0..3 {
                        out.push(equal_notes[i] as f64 + tuning_vector[i]);
                    }
                    out
                };
                if !within_epsilon_of(just_notes[2] - just_notes[1], ET_TO_JUST[4]) {
                    failures += format!("First inversion triad failed to tune perfectly: M3 was {}, but expected {}", 
                    just_notes[2] - just_notes[1], ET_TO_JUST[4]).as_str();
                }
                if !within_epsilon_of(just_notes[1] - just_notes[0], ET_TO_JUST[5]) {
                    failures += format!("First inversion triad failed to tune perfectly: p4 was {}, but expected {}", 
                    just_notes[1] - just_notes[0], ET_TO_JUST[4]).as_str();
                }
            }

            assert!(failures.is_empty(), "{failures}");
        
        }

        #[test]
        fn cluster_test() {
            //stress test: ensure that dense tone clusters don't have
            //too much drift from the original pitch centers.
            //Mostly of theoretical interest: no performer cares if their tone clusters are in tune!

            let DRIFT_TOLERANCE = 0.1; 

            //In practice, we find that drift tends to cap out around 7 cents.

            //Interestingly, drift seems to cap out at 4-note clusters: why could this be?
            //probably some cool linear algebra at play here :D

            let mut failures = String::new();

            //Whole step clusters
            for root in 11..80 {
                for cluster_size in 1..9 {
                    let equal_cluster = {
                        let mut out = vec![];
                        for i in 0..cluster_size {
                            out.push(root + 2 * i);
                        }
                        out
                    };
                    let tuning_vector = compute_tuning_vector_f64(&equal_cluster).unwrap();
                    for i in 0..tuning_vector.len() {
                        if f64::abs(tuning_vector[i]) >= DRIFT_TOLERANCE {
                            failures += format!("Too much drift in whole step cluster! Drift = {}, with {} notes in cluster. \n", tuning_vector[i], cluster_size).as_str();
                        }
                    }
                }    
            }

            //Half step clusters
            for root in 11..80 {
                for cluster_size in 1..9 {
                    let equal_cluster = {
                        let mut out = vec![];
                        for i in 0..cluster_size {
                            out.push(root + i);
                        }
                        out
                    };
                    let tuning_vector = compute_tuning_vector_f64(&equal_cluster).unwrap();
                    for i in 0..tuning_vector.len() {
                        if f64::abs(tuning_vector[i]) >= DRIFT_TOLERANCE {
                            failures += format!("Too much drift in half step cluster! Drift = {}, with {} notes in cluster. \n", tuning_vector[i], cluster_size).as_str();
                        }
                    }
                }    
            }

            assert!(failures.is_empty(), "{failures}");

            //Half step clusters

        }
    }
    
    mod auditory_tests {

        //uses simple normalized triangle waves to play chords.

        use super::*;
        mod wavmaker {
            use hound;

            use super::ET_TO_JUST;

            const CHANNELS: u16 = 1;
            const SAMPLE_RATE: u32 = 48000;
            const BIT_DEPTH: u16 = 16;
            const SPEC: hound::WavSpec = hound::WavSpec {
                channels: CHANNELS,
                sample_rate: SAMPLE_RATE,
                bits_per_sample: BIT_DEPTH,
                sample_format: hound::SampleFormat::Int
            };

            const PITCH_FREQ_CONVERSION_CONSTANT: f64 = 8.175798915643707;
            //pre-computed value of 440.0 / f64::powf(2.0, 69.0 / 12.0);

            fn pitch_to_freq(note: f64) -> f64 {
                return f64::powf(2.0, note / 12.0) * PITCH_FREQ_CONVERSION_CONSTANT;
            }
            fn make_triangle_wave(note: f64, duration_seconds: f64) -> Vec<f64> {

                let freq = pitch_to_freq(note);
                let total_samples = (duration_seconds * SAMPLE_RATE as f64) as usize;
                let period_in_samples = (SAMPLE_RATE as f64 / freq) as usize;
                let slope_in_samples = 4.0 / period_in_samples as f64;

                let mut out = Vec::with_capacity(total_samples);

                for i in 0..total_samples {
                    let phase = i % period_in_samples;
                    let current_sample = match phase {
                        p if p < period_in_samples / 4 => {
                            slope_in_samples * p as f64
                        },
                        p if p < (period_in_samples * 3) / 4 => {
                            2.0 - slope_in_samples * p as f64
                        },
                        p => {
                            -4.0 + slope_in_samples * p as f64
                        }
                    };
                    out.push(current_sample);
                }
                out
            }

            fn export_wav_from_vec_of_samples(filepath: &str, samples: Vec<f64>) {
                let mut wav_writer = hound::WavWriter::create(filepath, SPEC).unwrap();
                for sample in samples {
                    wav_writer.write_sample((sample * i16::MAX as f64) as i16).unwrap();
                }
                wav_writer.finalize().unwrap();
            }

            fn make_chord(notes: Vec<Vec<f64>>) -> Vec<f64> {
                //adds all, then normalizes.
                let mut max_len = 0;
                for note in  &notes {
                    if max_len < note.len() {
                        max_len = note.len();
                    }
                }

                let mut out = Vec::with_capacity(max_len);

                for note in &notes {
                    for i in 0..note.len() {
                        if i < out.len() {
                            out[i] += note[i];
                        } else {
                            out.push(note[i]);
                        }
                    }
                }

                let normalizer = notes.len() as f64;

                for sample in &mut out {
                    *sample = *sample / normalizer;
                }
                out
            }

            mod wav_creation_tests {
                use super::*;
                //Not needed to test the algorithm itself:
                //These were used early in testing to ensure the WAV production was accurate.
                //#[test]
                fn test_wav_construction() {
                    let a4_triangle = make_triangle_wave(69.0, 2.0);
                    export_wav_from_vec_of_samples("./test_A4.WAV", a4_triangle);
                    //user must manually listen to see if it's any good:
                    //in test suite so it won't compile in build versions where it's not needed.
                }

                //#[test]
                fn test_et_chord_construction() {
                    let duration = 2.0;
                    //cm7
                    let mut notes = vec![];
                    notes.push(make_triangle_wave(60.0, duration)); //C4
                    notes.push(make_triangle_wave(63.0, duration)); //Eb4
                    notes.push(make_triangle_wave(67.0, duration)); //G4
                    notes.push(make_triangle_wave(70.0,duration )); //Bb4
                    export_wav_from_vec_of_samples("./test_Cm7.WAV",
                        make_chord(notes)
                    );
                }

                //#[test]
                fn test_just_table() {
                    //This might be too subtle for technical people / non-musicians to hear, though. :(
                    //using the table, not the automatic tuner
                    let duration = 2.0;
                    //cm7
                    {
                        let mut notes = vec![];
                        notes.push(make_triangle_wave(60.0, duration)); //C4
                        notes.push(make_triangle_wave(60.0 + ET_TO_JUST[4], duration)); //E4
                        //notes.push(make_triangle_wave(60.0 + ET_TO_JUST[7], duration)); //G4
                        export_wav_from_vec_of_samples("./test_Cmaj_just_table.WAV",
                            make_chord(notes)
                        );
                    }
                    {
                        let mut notes = vec![];
                        notes.push(make_triangle_wave(60.0, duration)); //C4
                        notes.push(make_triangle_wave(64.0, duration)); //E4
                        //notes.push(make_triangle_wave(67.0, duration)); //G4
                        export_wav_from_vec_of_samples("./test_Cmaj_et.WAV",
                            make_chord(notes)
                        );
                    }
                }
            }
        }   
        
        #[test]
        fn stacked_perfect_interval_tests() {
            //p4 x5
            let max = 3;
            let mut equal_notes = vec![0];
            for i in 0..max {
                equal_notes.push(equal_notes[i] + 5);
            }
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("Stacked 4ths: notes: {:?}", equal_notes);
            println!("tuning: {:?}", tuning_vector);
    
            //p5 x5
            let mut equal_notes = vec![0];
            for i in 0..max {
                equal_notes.push(equal_notes[i] + 7);
            }
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("stacked 5ths: notes: {:?}", equal_notes);
            println!("tuning: {:?}", tuning_vector);
        }
    
        //non-rigorous-- sanity check to see how the algorithm affects more complicated chords
    
        fn extension_tests() {
            //dom7
            let equal_notes = vec![0, 4, 7, 10];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("dominant 7: notes: {:?}", equal_notes);
            println!("tuning: {:?}", tuning_vector);
            //lowers the fifth by a ton: not ideal!
    
            //m9
            let equal_notes = vec![0, 7, 10, 14, 15];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("minor 9: notes: {:?}", equal_notes);
            println!("tuning: {:?}", tuning_vector);
            //seems normal
    
            //M7#11
            let equal_notes = vec![0, 7, 11, 16, 18];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes).unwrap();
            println!("maj7#11: notes: {:?}", equal_notes);
            println!("tuning: {:?}", tuning_vector);
            //seems normal
        }
    
        fn asymmetry_tests() {
           //test on things like dominant 7th chords with 2 roots represented, to ensure it doesn't drag up the 7th too much.
           
        }
    }
    
}

fn main() {

}