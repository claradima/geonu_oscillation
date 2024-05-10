Most of the files here are not useful, I just uploaded them as backup/to keep track of what I did. There are three important files:

1) Earth models more bins - 2.ipynb  -> This explains how the Earth model was implemented in detail and it includes sources for all the data and formulas I used
   
2) Survival probability check corrected sigma IBD.ipynb  -> This is the final code I ran to get all the plots in the presentation; The code for the Earth model is copied
   from the other file; the formula for sigma_IBD is different and the implementation for the survival probability is also a bit different, I corrected some of the
   oscillation parameters, but the method is pretty much the same

3) FULL survival probability check corrected sigma IBD.ipynb  -> Mostly same as the other, but full three-flavour P_ee used instead of the approximation; only ran this with 1D grid
   size 100 so far to check if it's all ok

The Earth grid is not fine enough but my laptop couldn't really handle more; I might need to rerun it when I get my new laptop (should happen very soon, will update the
README file and the plots in the presentation when I get it and rerun everything
