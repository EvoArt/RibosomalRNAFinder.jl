# Registration

ProfileHMMs 0.1.0 is registered in General (PR #176394, merged 18 November 2026).
FindRNA's CI now installs dependencies from General without a development override.

1. Run `julia --project=. -e 'using Pkg; Pkg.test()'` and check CI against
   the registered ProfileHMMs release.
2. Comment `@JuliaRegistrator register` on the tested FindRNA commit.
   New packages wait three days in General.
3. AutoMerge will flag the name FindRNA as close to the registered FoldRNA
   (Damerau–Levenshtein distance 2). Reply on the General PR explaining the
   difference (rRNA gene finding vs RNA folding) and ask for a manual merge.
4. TagBot tags each release once its registry PR is merged.
