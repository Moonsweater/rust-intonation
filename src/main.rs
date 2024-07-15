
//uses https://gregstrohman.com/wp-content/uploads/2019/10/jic01-Interval-tuning.pdf as guideline for interval-based JI tuning.

static et_to_just: [f64; 12] = [
    0.0, // 0, unis.
    1.1173, // 1, second
    2.0391, // 2, != 12 - et_to_just[10]
    3.1564, // 3
    3.8631, // 4
    4.9804, // 5
    5.8251, // 6, != 12 - et_to_just[6]
    7.0196, // 7
    8.1369, // 8
    8.8436, // 9
    9.6883, // 10, != 12 - et_to_just[2]
    10.8827 // 11
]; //Maps integer intervals mod 12 to just intervals

fn compute_tuning_vector_f64(equal_notes: &Vec<i8>) -> Vec<f64> {

    let just_intervals = equal_notes_to_just_intervals(equal_notes);

    let n = equal_notes.len();
    let k = (n-1) * n / 2;

    //Presumes the notes vector is sorted.

    // Least squares: choose tuning vector that minimizes
    // ||A * (tuning + equal) - just||^2,
    // Where A takes a vector of pitches to a vector holding
    // the interval from pitch 1 to 2, then 1 to 3, ... then 2 to 3, ... then n to n.
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
    //tuning = equal - (A^T * A)^-1 * A^T * et_to_just(A * equal)

    //We can show that A^T * A is n-1 on the diagonal and -1 elsewhere.
    //Additionally, (A^T * A)^-1  is 2/(n+1) on the diagonal and 1/(n+1) elsewhere.
    //AND, (A^T * A)^-1 * A^T is simply A^T * 1/(n+1).

    let mut tuning: Vec<f64> = vec![];
    
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
        if start_note == n-2 && end_note == n-1 && j != k-1 {
            panic!();
        }
        if j == k-1 && !(start_note == n-2 && end_note == n-1) {
            panic!();
        }
    }

    for i in 0..n {
        tuning[i] = tuning[i] / ((n+1) as f64) + equal_notes[i] as f64;
    }

    tuning
  
}

fn equal_notes_to_just_intervals(equal_notes: &Vec<i8>) -> Vec<f64> {

    //assumes equal_notes is sorted

    let n = equal_notes.len();
    let mut start_note = 0;
    let mut end_note = 1;
    let mut i = 0;

    let mut just_intervals = vec![];

    loop {
        just_intervals[i] = et_to_just[((equal_notes[end_note] - equal_notes[start_note]) % 12) as usize];
        i += 1;
        start_note += 1;
        end_note += 1;
        if end_note >= n {
            start_note += 1;
            end_note = start_note + 1;
        }
    }
}

#[cfg(test)]
//this thing is important, somehow. Something something, rust attributes.
mod tests {
    use super::*; //use everything in the surrounding module environment.
    fn unison_tests() {
        //ensure tuning_vector for single notes is zero
        for i in 0i8..127 {
            let equal_notes = vec![i];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            for j in 0..tuning_vector.len() {
                assert_eq!(0.0, tuning_vector[j])
            }
        }
        //ensure tuning vector for octaves is zero
        for i in 0i8..127 {
            let equal_notes = vec![i, (i + 12) % 120, (i + 12 * 4) % 120];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            for j in 0..tuning_vector.len() {
                assert_eq!(0.0, tuning_vector[j])
            }
        }
    }
    fn interval_tests(){
        //ensure tuning vector for dyads is as given by the table.
        let epsilon = 0.0001;
        for root in 0i8..100 {
            for top in 0i8..27 {
                let equal_notes = vec![root, top];
                let tuning_vector = compute_tuning_vector_f64(&equal_notes);
                for j in 0..tuning_vector.len() {
                    assert!((tuning_vector[j] - et_to_just[(top - root) as usize]) < epsilon)
                }
            }
        }
    }
    fn triad_tests() {
        //probably needs us to listen with our ears at this point :pensive:
        //major triad:
        //in 1st:
        for root in 0i8..100 {
            let equal_notes = vec![root, root + 4, root + 7];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            //play audio
        }
        //in 2nd:
        for root in 0i8..100 {
            let equal_notes = vec![root, root + 3, root + 8];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            //play audio
        }
        //in 3rd:
        for root in 0i8..100 {
            let equal_notes = vec![root, root + 5, root + 9];
            let tuning_vector = compute_tuning_vector_f64(&equal_notes);
            //play audio
        }
        
    }
    fn asymmetry_tests() {
        //test on things like dominant 7th chords with 2 roots represented, to ensure it doesn't drag up the 7th too much.
        
    }
}

fn main() {



}
