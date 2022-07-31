/*
SDUST complexity-filtering algorithm

Morgulis A, Gertz EM, Schäffer AA, Agarwala R. A fast and symmetric DUST implementation to mask
low-complexity DNA sequences. J Comput Biol. 2006 Jun;13(5):1028-40.
doi: 10.1089/cmb.2006.13.1028. PMID: 16796549.

This implementation is extracted from minimap2 by Heng Li:
https://github.com/lh3/minimap2/blob/master/sdust.h
the MIT license of which follows.

--

The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef SDUST_H
#define SDUST_H

#define SDUST_VERSION "0.1-r2"

#include <stdint.h>

struct sdust_buf_s;
typedef struct sdust_buf_s sdust_buf_t;

uint64_t *sdust(void *km, const uint8_t *seq, int l_seq, int T, int W, int *n);

#endif
