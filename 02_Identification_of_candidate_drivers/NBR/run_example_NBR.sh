ls Regions/*.bed | xargs -i Rscript find_noncoding_drivers_precomp_correctoverlaps.R mutations.txt {}
