set makeprg=[[\ -f\ Makefile\ ]]\ &&\ make\ \\\|\\\|\ make\ -C\ build/

set grepprg=grep\ -nrI\ --exclude-dir={.git,docs,build,libraries}\ --exclude=*tags

" recreate tags file with F5
map <F5> :GenerateVisarTags<CR>

set tags+=visar.tags
set tags+=opencv.tags

function GenerateVisarTags()
  !ctags -R --sort=yes --langmap=C++:+.cu --c++-kinds=+p --fields=+iaS --extra=+qf -f ./visar.tags .
endfunction
command GenerateVisarTags execute GenerateVisarTags()

function GenerateOpenCVTags()
  !ctags -R --sort=yes --langmap=C++:+.cu --c++-kinds=+p --fields=+iaS --extra=+qf --exclude=build --exclude=.git --exclude=cmake -f ./opencv.tags ~/software/opencv-2.4.13/
endfunction
command GenerateOpenCVTags execute GenerateOpenCVTags()
