#!/bin/sh

# Prints the quoted version information as a c(xx)flag directly.
# It can be used in a makefile via `<script-name>`.
echo -DVERSION=\"$(git describe --dirty)\"
