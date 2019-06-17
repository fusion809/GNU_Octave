function push(x)
  system "git add --all"
  system 'git commit -m "x"'
  system "git push origin $(git-branch)"
endfunction
