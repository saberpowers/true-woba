True wOBA: Estimation of true talent level for batters
======================================================

This repository contains the paper and slides for the presentation by Scott
Powers and Eli Shayer, titled ``True wOBA: Estimation of true talent level for
batters''  at the 2016 SABR Analytics Conference, as well as code for producing
the results in the paper.

In order to reproduce the results, you will need to use the function defined in
code/parse.retrosheet2.pbp.r (adapted from Marchi and Albert, 2014)
to automatically download play-by-play data from Retrosheet. To do so, simply
call

~~~R
parse.retrosheet2.pbp(2015)
~~~

in R, and then move download.folder/unzipped/all2015.csv to data/all2015.csv.
Then you can use code/regression.R to reproduce the results in Section 1 of the
paper and code/validation.r to reproduce the result in Sections 2 and 3.
