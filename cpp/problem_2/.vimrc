set number
set encoding=UTF-8
set tabstop=4
set softtabstop=4
set shiftwidth=4
set smarttab
set expandtab
set autoindent
set smartindent
set colorcolumn=100
highlight ColorColumn ctermbg=red
syntax on
filetype on
set cursorline
set showmatch
set title
colorscheme peachpuff
set mouse=a

call plug#begin()
Plug 'vim-airline/vim-airline'
Plug 'vim-airline/vim-airline-themes'
call plug#end()

let g:airline#extensions#tabline#enabled = 1
let g:airline#extensions#tabline#left_sep = ' '
let g:airline#extensions#tabline#left_alt_sep = '|'
let g:airline#extensions#tabline#formatter = 'unique_tail'
let g:arline_powerline_fonts = 1
let g:airline_theme='badwolf'
let g:airline#extensions#filetype#enabled = 1
let g:airline#extensions#modify#enabled = 1
let g:airline#extensions#tabline#enabled = 1
let g:airline#extensions#tabline#left_sections = ['branch']
