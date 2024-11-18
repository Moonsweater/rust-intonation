# Just Intonation Tuner
Library to automatically tune polyphonic signals to a best approximation of just intonation.
## Motivation
In ensemble music, performers naturally adjust their pitch to achieve mathematically "pure" intervals - what music theory calls just intonation. While the method for doing so is straightforward and well-documented when playing only two pitches at once, playing three or more simulatneous pitches complicates the matter, as sources only define just tuning between pairs of notes. This project formalizes and implements a method of achieving just intonation for 3 or more pitches at once, to enable software synthesizers to replicate human performers in this respect. 
## Formalism
This project approaches the task at hand as an optimization problem, seeking to minimize the weighted L^2 distance between each interval formed by a collection of simultaneous pitches and their just intonation equivalents.

Let `equal` be a vector of pitches (e.g, C2, A3, F#6, etc.), represented as integers, and sorted from lowest to highest. Let `A` be matrix. It can be shown that `A` is the matrix

```
    [-1, 1, 0, 0, 0, ..., 0, 0]
    [-1, 0, 1, 0, 0, ..., 0, 0]
    ...
    [0, -1, 1, 0, 0, ..., 0, 0]
    [0, -1, 0, 1, 0, ..., 0, 0]
    ...
    [0, 0, 0, 0, 0, ..., -1, 1]
```
.
It can be shown that `(A * x) mod 12` will map a vector `x` of pitches to a vector of every interval class formed by pairs of pitches in `x`.

We now seek to pick a vector of pitches `tuned` that minimizes the least-squares expression:
 
`||W((A * tuned) - just((A * equal) mod 12)||^2`,

Where `just` is the function that maps interval classes to their just-intonation counterparts, each measured in semitones.

`W`, conversely, is a diagonal matrix encoding the relative importance that each interval sound perfectly in-tune, musically speaking. `W` is determined at the time of computation by the exact intervals contained in `(A * equal) mod 12`, and in general, rates consonant intervals, like octaves or major thirds, more highly than dissonant intervals, like tritones or minor seconds, in keeping with a typical performer's way of thinking about just tuning.

