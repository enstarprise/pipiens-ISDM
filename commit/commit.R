hist(rnorm(1000,0,0.5))
# https://rfortherestofus.com/2021/02/how-to-use-git-github-with-r

library(usethis)
usethis::create_from_github(
  "https://github.com/enstarprise/pipiens-ISDM.git",
  destdir = "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/AIM 1/pipiens-ISDM/"
)


usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
# git add -A
# git commit -m "Ready for fresh GitHub push"
usethis::use_github()
# git remote -v  # Should show your NEW GitHub URL
# git push -u origin main  # Should work without 408 errors
# 
# git remote add origin https://github.com/enstarprise/pipiens-ISDM.git
# git branch -M main
# git push -u origin main