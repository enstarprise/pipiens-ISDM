# 1. Remove the corrupted .gitignore
rm .gitignore

# 2. Create a new, clean .gitignore with proper formatting
cat > .gitignore << 'EOF'
# Ignore ALL files by default
*
  
  # Then WHITELIST specific file types
  !*.R
!*.r
!*.csv
!.gitignore
!README*
  !*.md
!*.Rmd
!*.qmd

# Whitelist directory structure (keeps empty directories)
!*/
  
  # Optional: Whitelist specific directories
  !src/
  !R/
  !data/
  !output/
  
  # Re-ignore specific file types within whitelisted directories
  output/*.pdf
output/*.png
data/raw/
  *.pdf
*.png
*.jpg
*.zip
*.rds
*.rda
*.RData
EOF



using the R terminal
# 1. Remove everything from Git tracking (keeps files locally)
git rm -r --cached .

# 2. Add your clean .gitignore
git add .gitignore

# 3. Stage only files that match your whitelist
git add -A

# 4. Check what will be committed (should show only .R, .r, .csv files)
git status

# 5. Commit the cleanup
git commit -m "Clean repo: keep only .R, .r, and .csv files"


cat .gitignore

# 6. Push the much smaller repository
git push origin main



