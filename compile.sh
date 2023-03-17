#!/bin/bash
#
# Compilation script.

if [ ! -d builddir ]; then
	mkdir builddir
	meson setup builddir
fi

cd builddir
meson compile

cd ..
for file in builddir/bayes.cpython-*.so; do
	cp "$file" .
done