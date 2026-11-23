# Registration

The package was first submitted to General as FindRNA ([General#177276](https://github.com/JuliaRegistries/General/pull/177276)).
AutoMerge blocked that name as too close to the registered FoldRNA
(Damerau–Levenshtein distance 2), so it was renamed RibosomalRNAFinder, keeping the
same UUID. No registered name is within distance 3 of RibosomalRNAFinder.

1. Run `julia --project=. -e 'using Pkg; Pkg.test()'` and check CI.
2. Comment `@JuliaRegistrator register` on the tested commit. New packages wait
   three days in General.
3. TagBot tags each release once its registry PR is merged.
