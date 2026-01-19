#!/bin/sh

set -xe

# Linux used for mainly debugging

if [ "$1" = "all" ]; then
    CC=gcc
    CFLAGS="-g -O0 -Wall -fsanitize=address -Wextra -fdebug-prefix-map=$(pwd)=. `pkg-config --cflags raylib`"
    LIBS="`pkg-config --libs raylib` -lglfw -ldl -lpthread -lm"
    OUT=spectrogram
    $CC $CFLAGS spectrogram.c -o $OUT $LIBS
    CC=x86_64-w64-mingw32-gcc
    CFLAGS="-Wall -Wextra -I$HOME/dev/raylib/src"
    LIBS="-L$HOME/dev/raylib/src -lraylib -lopengl32 -lgdi32 -lwinmm"
    OUT=spectrogram.exe
    $CC $CFLAGS spectrogram.c -o $OUT $LIBS
else
    if [ "$1" = "linux" ]; then
        CC=gcc
        CFLAGS="-g -O0 -Wall -fsanitize=address -Wextra -fdebug-prefix-map=$(pwd)=. `pkg-config --cflags raylib`"
        LIBS="`pkg-config --libs raylib` -lglfw -ldl -lpthread -lm"
        OUT=spectrogram
    else
        CC=x86_64-w64-mingw32-gcc
        CFLAGS="-Wall -Wextra -I$HOME/dev/raylib/src"
        LIBS="-L$HOME/dev/raylib/src -lraylib -lopengl32 -lgdi32 -lwinmm"
        OUT=spectrogram.exe
    fi 
        $CC $CFLAGS spectrogram.c -o $OUT $LIBS
fi