
//uses https://gregstrohman.com/wp-content/uploads/2019/10/jic01-Interval-tuning.pdf as guideline for interval-based JI tuning.

static ET_TO_JUST: [f64; 12] = [
    0.0, // 0, unis.
    1.1173, // 1, second
    2.0391, // 2, != 12 - ET_TO_JUST[10]
    3.1564, // 3
    3.8631, // 4
    4.9804, // 5
    5.8251, // 6, != 12 - ET_TO_JUST[6]
    7.0196, // 7
    8.1369, // 8
    8.8436, // 9
    9.6883, // 10, != 12 - ET_TO_JUST[2]
    10.8827 // 11
]; //Maps integer intervals mod 12 to just intervals

fn compute_tuning_vector_f64(equal_notes: &Vec<i8>) -> Vec<f64> {

    let n = equal_notes.len();
    // let equal_intervals = {let mut a = vec![];
    //     for i in 0..n {
    //         a.push(equal_notes[i] - equal_notes[0]);
    //     }
    //     a
    // };
    if n == 0 {return vec![];}
    if n == 1 {return vec![0.0];}
    let k = ((n-1) * n ) / 2;

    let just_intervals = equal_notes_to_just_intervals(equal_notes);

    //Presumes the notes vector is sorted.

    // Least squares: choose tuning vector that minimizes
    // ||A * (tuning + equal) - just||^2,
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

    //This makes our least squares problem a quest to solve:
    //A * tuning = A(A^T * A)^-1 * A^T * (A * equal - just), i.e.,
    //tuning = equal - (A^T * A)^-1 * A^T * just, i.e.,
    //tuning = equal - (A^T * A)^-1 * A^T * ET_TO_JUST(A * equal)

    //We can show that A^T * A is n-1 on the diagonal and -1 elsewhere.
    //Note, however, that A^T * A is singular. This messes with the least squares.
    //To work around this, we use the pseudoinverse of A^T * A. This gives us:
    //A(tuning + equal) = A * (A^T * A)^(+) * A^T * just
    // ==> tuning = A^(+) * just - equal
    //The pseudoinverse of A can be shown to be A^T / n, due to all its nonzero singular values being n.

    //Since the nullspace of A^T * A is span([1; 1; ... 1]), the best answer will be A^(T) * just - equal, plus some scalar multiple of [1; 1; ... 1].
    //We choose a linear shift that makes the bass note perfectly in tune with equal temperment.
    //(More sophisticated implementations would identify the *root*, and make *that* in tune.)

    let mut tuning: Vec<f64> = vec![0.0; n];
    
    //will rustify later with scoping and stuff.

    let mut start_note: usize = 0;
    let mut end_note: usize = 1;

    for j in 0..k {
        tuning[start_note] -= just_intervals[j];
        tuning[end_note] += just_intervals[j];
        end_note += 1;
        if end_note >= n {
            start_note += 1;
            end_note = start_note + 1;
        }
        if start_note >= n || end_note >= n{
            println!("Ending at start_note = {}, j = {}.", start_note, j);
            break;
        }
        // if start_note == n-2 && end_note == n-1 && j != k-1 {
        //     panic!();
        // }
        // if j == k-1 && !(start_note == n-2 && end_note == n-1) {
        //     panic!();
        // }
    }

    for i in 0..n {
        tuning[i] = (tuning[i] / (n as f64)) - equal_notes[i] as f64;
    }
    //linear offset by multiple of [1; ... 1], so that tuning[0] = 0.
    let offset = tuning[0];
    for i in 0..n {
        tuning[i] -= offset;
    }

    tuning
  
}

fn equal_notes_to_just_intervals(equal_notes: &Vec<i8>) -> Vec<f64> {

    //assumes equal_notes is sorted and has len >= 2

    let n = equal_notes.len();
    let mut start_note = 0;
    let mut end_note = 1;

    let mut just_intervals = vec![];

    loop {
        just_intervals.push(
            {let interval = 
            {let octaves = (((equal_notes[end_note] - equal_notes[start_note]) / 12) * 12) as f64;
            println!("{} octaves", octaves / 12.0);
            octaves}
            + ET_TO_JUST[((equal_notes[end_note] - equal_notes[start_note]) % 12) as usize];
            println!("interval equals {}", interval);
            interval
            }
        );
        
        end_note += 1;
        if end_note >= n {
            start_note += 1;
            end_note = start_note + 1;
        }
        if start_note >= n || end_note >= n{
            break;
        }
    }

    just_intervals

}

#[cfg(test)]
mod tests {
    use super::*;
    const TEST_LO: i8 = 36; //C2
    const TEST_HI: i8 = 100; //E7
    fn jury_rigged_random(a: i8) -> i8 {
        if rand::random() && a < 4 {
            jury_rigged_random(a + 1)
        } else {
            a
        }
    }
    #[test]
    fn unison_tests() {
        //ensure tuning_vector for single notes is EMPTY (not just zero)
        println!("Unison:");
        for i in TEST_LO..TEST_HI {
            let equal_notes = vec![i];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            for j in 0..tuning_vector.len() {
                println!("tuning_vector[{}] = {}", j, tuning_vector[j]);
                assert_eq!(0.0, tuning_vector[j]);
            }
        }
        //ensure tuning vector for octaves is zero
        println!("Octaves:");
        for i in TEST_LO..(TEST_HI / 2) {
            let r1 = jury_rigged_random(0);
            let r2 = jury_rigged_random(1);
            let equal_notes = 
                {let mut a = vec![i, (i + 12 * r1), (i + 12 * r2)];
                a.sort_unstable(); a};
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            for j in 0..tuning_vector.len() {
                println!("equal_notes = [{}, {}, {}]", equal_notes[0], equal_notes[1], equal_notes[2]);
                println!("tuning_vector = [{}, {}, {}]", tuning_vector[0], tuning_vector[1], tuning_vector[2]);
                assert_eq!(0.0, tuning_vector[j])
            }
        }
    }

    #[test]
    fn interval_tests(){
        //ensure tuning vector for dyads is as given by the table.
        println!("Intervals:");
        let epsilon = 0.0001;
        for root in TEST_LO..TEST_HI {
            for interval in 0i8..27 {
                let top = root + interval;
                let equal_notes = vec![root, top];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes);
                for j in 0..tuning_vector.len() {
                    assert!((tuning_vector[j] - ET_TO_JUST[(interval % 12) as usize]) < epsilon)
                }
            }
        }
    }
    
    //fn triad_tests() {
    //    //probably needs us to listen with our ears at this point :pensive:
    //    //major triad:
    //    //in 1st:
    //    for root in TEST_LO..TEST_HI {
    //        let equal_notes = vec![root, root + 4, root + 7];
    //        let tuning_vector = compute_tuning_vector_f64(&equal_notes);
    //        //play audio
    //    }
    //    //in 2nd:
    //    for root in TEST_LO..TEST_HI {
    //        let equal_notes = vec![root, root + 3, root + 8];
    //        let tuning_vector = compute_tuning_vector_f64(&equal_notes);
    //        //play audio
    //    }
    //    //in 3rd:
    //    for root in TEST_LO..TEST_HI {
    //        let equal_notes = vec![root, root + 5, root + 9];
    //        let tuning_vector = compute_tuning_vector_f64(&equal_notes);
    //        //play audio
    //    }
    //    
    //}
    //fn asymmetry_tests() {
    //    //test on things like dominant 7th chords with 2 roots represented, to ensure it doesn't drag up the 7th too much.
    //    
    //}
    
    
}

fn main() {

}
