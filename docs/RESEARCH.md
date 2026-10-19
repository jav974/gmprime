# Research log

How gmprime's coding toolkit and gmz came to be: what was tried, what it gave, and what was kept or set aside.
Every size below is a real compressed message, decompressed and compared with the original. Test files come
from the [Canterbury corpus](https://corpus.canterbury.ac.nz/), plus `source.txt` (125,000 bytes of this
project's sources).

## Reference points

| File | Size | gzip -9 | bzip2 -9 | xz -9e | zstd -22 | brotli -11 |
|---|---|---|---|---|---|---|
| alice29.txt | 152,089 | 54,194 | 43,202 | 48,536 | 49,215 | 46,487 |
| asyoulik.txt | 125,179 | 48,803 | 39,569 | 44,600 | 45,141 | 42,713 |
| lcet10.txt | 426,754 | 144,462 | 107,706 | 119,496 | 121,340 | 113,416 |
| plrabn12.txt | 481,861 | 194,357 | 145,577 | 165,464 | 167,442 | 163,267 |
| source.txt | 125,000 | 26,938 | 23,555 | 23,824 | 24,326 | 22,846 |

On English prose, bzip2 (Burrows-Wheeler) is the best of the common tools. Research context-mixing compressors
(the PAQ family) reach about 1.7 bits per character on alice29.txt.

## Compression

### Measuring first

- Order-k entropy of alice29.txt falls from 4.91 bits/byte (order 0) to 1.93 (order 2) and 0.62 (order 5), but an
  adaptive model actually coding the data does best at order 2 (52 KB) and worse above: contexts get too rare to
  learn. Entropy estimates are optimistic; `adaptive_code_length` is the honest figure.
- 95% of source.txt is covered by 8,498 repeats of 4 bytes or more. On 40 copies of it, 7 more repeats cover
  99.9%: what xz exploits (24.6 KB) and no context model sees (988 KB at best).

### From one order to context mixing (alice29.txt)

| Step | Bytes |
|---|---|
| one order (2), adaptive | 52,818 |
| orders 0 to 6, mixed by a logistic mixer | 41,559 |
| mixer weights selected by the partial byte | 40,115 |
| + match model | 39,777 |
| + probability map (SSE) and word model | 38,094 |
| + adaptive order: weights selected by the highest order already seen | 37,615 |
| + learning rates tuned | 37,426 |
| + indirect modeling (bit histories mapped to probabilities) | 37,073 |

Dead ends: two counters per context at different speeds (worse: more inputs for the mixer), a second probability
map on order 2 (no gain).

### Other views of the data

| View | alice29.txt | Why |
|---|---|---|
| the file as a decimal number | 152,745 (incompressible) | base conversion mixes all bits: the value of a file has no structure, its bit positions do |
| file XOR its prediction (55% right) | 59,455 | a hard prediction throws away its uncertainty; models must give probabilities |
| bit planes, plane by plane, past only | 59,146 | |
| bit planes with the higher planes of the following bytes known | 45,755 | seeing the future helps a lot, but byte order still wins |
| skeleton (top 3 planes) first, then bytes with the future skeleton | 41,827 | coding the skeleton without the past bytes costs what the future saves |
| plain byte order, same contexts | 38,394 | |

The order in which bits are coded only moves the difficulty (chain rule); left to right is the order models learn
best for text.

### Parallel research directions

Five directions were explored in parallel against the same baseline (307,403 bytes on the five texts):

| Direction | Best | Kept |
|---|---|---|
| Global analysis sent as side information (repeat hints, priming, per-block parameters) | -0.5% at best, matched by a causal mixture of models | no: describing statistics costs as much as learning them, and mixing learns almost for free |
| Text: column contexts, capitals transform, mixer by character class, word dictionary | -1.76% | yes, but the dictionary |
| Views and masks: 2D line structure, what followed the last occurrence, periodicity, BWT, Gray code, bit reversal | -1.1 to -2.0% | line structure and last occurrences; every transform lost 2 to 8% |
| Model strength: slot checksums, deterministic inputs, more mixers, chained maps, rate schedule, runs | -1.86% | yes, but runs (they help code only) |
| Folding the bit field (2D, 3D, spiral, Hilbert) | fax image -50%, spreadsheet -54% | yes, for data that folds; spiral and Hilbert orders lost 9 to 39% |

Folding is the striking one: `popcount((N >> L) ^ N)` over the whole file as one number finds the row width of a
fax image (1728 bits) and the record size of a spreadsheet (104 bits) by itself, and templates across the fold
halve what xz achieves on them. Per-file parameter search gained only 0.06 to 0.32%, but showed that slower counters
suit every file.

### gmz

Combining all that was kept, with `analyze_folds` choosing between the text and the folded model:

| File | bzip2 -9 | xz -9e | gmz |
|---|---|---|---|
| alice29.txt | 43,202 | 48,536 | 35,562 (1.87 bits/byte) |
| lcet10.txt | 107,706 | 119,496 | 86,054 |
| ptt5 | 49,759 | 39,868 | 21,860 |
| kennedy.xls | 130,280 | 51,876 | 12,968 |

The biggest contributors on text, measured by removing each: word contexts (6,018 bytes), the adaptive-order
mixer (2,782), chained maps (2,407), the learning-rate schedule (1,524).

### Enumerative coding

The data as the result of a count rather than of a prediction: knowing the length n and the k bits set, the
field is one of C(n, k) arrangements, coded by its rank (`combination_rank`). Exact rank in one pass, the current
binomial updated by a small multiplication and an exact division per bit: quadratic time, 0.2 s for 10^6 bits
at density 0.01 but 18 s at 0.5. Blocks of 4,096 bits ranked separately and packed in mixed radix
(`rank_blocks`) take 0.03 to 0.08 s, 0.03% to 0.85% above log2 C(n, k).

Masks of 10^6 bits, sizes in bits (k included):

| Mask | k | log2 C(n, k) | Rank | Adaptive order 0 | Adaptive order 12 |
|---|---|---|---|---|---|
| random | 100 | 1,488 | 1,487 | 1,640 | 1,672 |
| random | 10,000 | 80,805 | 80,804 | 81,176 | 81,544 |
| clustered | 1,000 | 11,421 | 11,412 | 9,264 | **5,632** |
| clustered | 100,000 | 469,006 | 469,005 | 465,000 | **369,416** |

On real files, in bytes: bit planes ranked in blocks give alice29.txt 110,820, ptt5 158,771; the arrangement of
the bytes given their counts (`multiset_rank`) alice29.txt 86,823, kennedy.xls 452,542; gmz 35,563 and 12,970.

The rank is optimal for given rules and pays nothing to learn, 0.4% to 10% better than adaptive coding on random
positions, most on small streams. But it sees only the counts: any structure beyond them goes to context models,
by far. Uses: small sparse side streams (exception positions, flags), sets and permutations (a dictionary), and
it is the one tool so far where big-number arithmetic is the algorithm itself.

### Broad evaluation

Calgary (18 files), Silesia (12 files, 212 MB) and the first 10 MB of enwik8, every gmz size a verified round
trip (gmz before the integer arithmetic and the machine code models), against paq8px v217 at -8 and zpaq at
-m5. Bits per byte:

| Corpus | brotli -q 11 | xz -9e | zpaq -m5 | paq8px -8 | gmz |
|---|---|---|---|---|---|
| Calgary | 2.108 | 2.174 | 1.772 | 1.480 | 1.680 |
| enwik8, 10 MB | 2.164 | 2.176 | 1.673 | 1.420 | 1.659 |
| Silesia | 1.871 | 1.829 | 1.476 | (too slow) | 1.539 |

- gmz is smaller than every standard tool on every file, by 6% to 45%.
- Against zpaq it wins on text (nci -7.3%, reymont -3.3%, samba -2.5%) and loses on ooffice (+21.8%), sao
  (+12.8%), mozilla (+6.1%) and osdb (+4.5%).
- paq8px stays 13% to 17% ahead on text and 33% on xml, at 11 to 19 KB/s against gmz's 90 to 135 KB/s for a
  round trip. On Silesia it took over 20 minutes for some files, so only xml and x-ray were completed (x-ray:
  gmz 4.4% larger).

Where gmz lost, the model did not fit the data:
- **Fold detection is too strict.** Forced fold mode gives sao -11.4% (records of 28 bytes, now below zpaq),
  mr -8.9% (it reaches the 8192-bit limit), geo and obj2 -4.2%; fold mode is almost never worse than text.
- **x86 code had no model.** The E8/E9 transform alone gains 15% on ooffice.
- **Silesia's mozilla is not x86** but Alpha (Tru64) binaries: no branch transform applies.

### Machine code

Gains on a 2 MB slice of ooffice (Windows DLLs, x86) and of x86-64 code, each step added to the previous ones:

| Step | ooffice | x86-64 |
|---|---|---|
| E8/E9 transform (reversible: encoding from the end, decoding from the start) | -14.9% | -27.8% |
| instruction contexts: field and index in the instruction, opcode, ModRM, previous opcodes | -9.2% | -12.6% |
| 3 previous opcodes; previous opcode with its ModRM | -1.7% | |
| indirect: the byte seen last time at the same place in an instruction | -1.4% | -2% |
| 5 sparse contexts (bytes 2 and 3 back alone, high nibbles, bytes 4 and 8 back, 32-bit alignment) | -1.5% | |
| mixers selected by the instruction state and by the opcode, a map on the instruction state | -0.6% | |

The text contexts then gain 0.03% for twice the time on x86 code: they are left out there. Also measured:
- transforming conditional jumps too loses 1.2%, as zpaq and xz found;
- the minimum match length (4, 6, 9) and the counter limits (127, 1023) make no difference.

The full ooffice drops from 2,151,593 to 1,542,295 bytes (-28.3%, zpaq 1,766,663), in 1.4 times the time.

Sparse contexts also help binary data that is not code: sum -8.9%, geo -8.3%, obj2 -8.1% (68k code, still 29%
behind paq8px), obj1 -5.6%, osdb -4.2%, sao -2.0%, mr -1.7%, x-ray -1.0%, for 30% to 60% more time.

Detection, never wrong on text:
- **binary:** 10% of control bytes. At 1%, samba (mostly source text) paid 35% more time for 0.2%.
- **x86:** 50 near CALL or JMP per 10,000 bytes (x86 code 280 to 313, Alpha 11, kennedy.xls 4, text 0).
- **64-bit:** frequent REX.W prefixes.

Forcing a wrong guess costs alice29.txt +6.2% as x86 and +0.5% as binary. A threshold still splits poorly
(trans at 9% gains 1%, a samba slice at 14.5% gains 0.3%): comparing on the sample, as for capitals, would.

### Records, lattices and images

The broad evaluation showed fold detection missing records and images. Measured first:
- `analyze_folds` finds 224 bits on sao (the 28-byte records), 32 bits on geo and 8192 bits (its limit) on mr,
  but none passes the strength threshold; x-ray's rows (30,400 bits) are beyond the limit.
- Byte autocorrelation shows sharp peaks: sao 28, geo 4, kennedy.xls 13, x-ray 2 and 3,800, mr 2 and 1,024.
- The entropy of each byte column of sao, given the byte a record above: 11 columns stay at 8 bits, the
  mantissas of its doubles.

**Choosing by comparison.** The analysis suggests ways (text, capitals, fold, bitmap, records, lattice, signal,
16-bit image), each compresses the first 64 KB, and the smallest message decides. A folded way must beat text by
1%, since it costs about twice the time. Only binary data or a strong fold goes this way: on text, trying the
structure plans tripled the time of cp.html for nothing.

**Record contexts.** The column of the byte in its record, with the byte above, two above, the previous byte, the
northeast and northwest ones, and the extrapolation 2N - NN; a mixer selected by the column and bit position.
sao -5.2% and geo -4.4% alone.

**The decimal lattice.** sao's doubles are integers in disguise: right ascension in milliseconds of time,
declination in hundredths of an arcsecond, proper motions in 0.0001 s and 0.001". Each field becomes k and a
residual of a few units in the last place (no exact formula for the residual was found: at best 14,000 of
71,000 values). k is written high byte first in the field's own bytes, so that its constant high bytes cost
nothing, with a byte of flags per record for values off the lattice. Decoding computes `k * step` with integers
only, exactly as IEEE multiplication rounds; 2 million random products checked against the processor. The same
tool covers geo's IBM floats (k / 64).

**16-bit images.** Byte order and row width are detected, pixels are coded high byte first, predicted from W, N,
their median, the gradient, 2N - NN and two integer LMS predictors; each bit's contexts hold where each
prediction stands against the middle of the values the known bits still allow, on a logarithmic scale, so that
statistics are shared across intensities.

Full Silesia files:

| File | gmz before | gmz | zpaq -m5 | paq8px -8 |
|---|---|---|---|---|
| sao | 4,399,334 | **2,047,265** | 3,899,367 | 3,727,085 |
| mr | 2,192,035 | **1,708,413** | 2,181,419 | |
| x-ray | 3,651,981 | **3,472,470** | 3,669,812 | 3,530,564 |

On sao, paq8px takes 22 minutes and 2.4 GB for a message 82% larger than gmz's, made in 54 s.

geo goes from 45,525 to 41,082 bytes (paq8px 41,988), sum -5.9%, obj2 -2.1%; Calgary and Canterbury together
-0.64%, no file worse by more than a byte. osdb, nci and kennedy.xls do not move. Compressing a structured file
takes 1.5 to 2 times longer, and small binary files more: every way compresses the whole file when it is under
64 KB.

Also tried: geo as a 2D image of traces (worse), more LMS weights (no change), and the per-byte period context,
which predicts each byte value back at its own last distance: -3.7% on sao and -1.9% on geo before the lattice,
nothing once it is in.

Integrating it found two bugs that macOS release builds hid:
- the image model wrote 12 contexts into a table of 11; the twelfth was never read, but the write landed
  outside the table;
- an 8-byte field without residual bits, in offset binary, shifted a signed 1 into the sign bit.

AddressSanitizer on Linux caught the first; the unit tests the second.

### The lattice on real floating point files

Is sao a lucky case? 23 public files that store values in binary floating point, each round trip verified:
shapefiles and a GeoPackage (TIGER roads, Natural Earth), FITS and raw tables of the Hipparcos catalog, NetCDF,
HDF5, SQLite, Parquet and Feather of a weather buoy, an Argo float profile, meshes, and negative controls (float
audio, satellite images, simulation output).

The lattice wins on decimal values of high cardinality, where almost every value is new: their mantissas' low
bytes look random to any model, and only the lattice makes them computed. It gains nothing on sensor readings of
low cardinality (a buoy: 100 to 500 distinct values per column), whatever the container: column and order
contexts already learn each value's bytes, and paq8px stays 1% to 55% ahead there. It is never chosen wrongly
on the negative controls.

Why the first version missed some:
- a header before the records (4 KB in the NetCDF, 11 KB in the FITS) sent the step search astray;
- values off the lattice mixed in (computed points among the places of Natural Earth);
- 8 fields at most (the buoy's records have 12 floats);
- records of varying length (polylines, GeoPackage and SQLite pages): no fixed width at all;
- pure streams of doubles taken for text (under 10% of control bytes), so never analyzed.

The fixes: decimal steps (powers of 10) tried first, robust to a quarter of the values off the lattice; fields
searched on records away from the start, up to 24; and the scan transform, which finds decimal doubles and floats
anywhere and codes each as its k, or its difference with the value one or two before (x and y alternate), as a
token or in place. Both the scan and the record plans are compared on the sample, the scan with a margin of 2%
(the sample may not show how scanned columns spread over the file).

| File | zpaq -m5 | gmz before | gmz | paq8px -8 |
|---|---|---|---|---|
| Natural Earth places (.shp) | 101,022 | 73,856 | **56,111** | 84,786 |
| TIGER roads (.shp) | 253,907 | 200,550 | **161,136** | (128 KB slice: gmz 33% below) |
| TIGER roads (.gpkg) | 280,582 | 235,236 | **170,105** | (slice: gmz 37% below) |
| TIGER coordinates as raw doubles | 195,776 | 203,108 | **115,645** | (slice: gmz 43% below) |
| Hipparcos (FITS table) | 307,166 | 238,122 | **223,704** | (slice: gmz 1% above) |
| buoy (SQLite) | 178,296 | 166,118 | 158,998 | (slice: gmz 55% above) |

Not done: choosing the step value by value for columns that mix precisions (USGS coordinates: 4 decimals and
1/60000 of a degree), and fields across the boundary of a record.

### Numbers in text

Logs, CSV, JSON and tables spend most of gmz's bits on digits: 61% to 98% of the coded bits for 15% to 74% of the
bytes. Every general compressor, paq8px included, codes them digit by digit as text.

**The transform.** Each number becomes a marker byte (absent from the file), its layout (sign, digits before and
after the point, leading zeros included) and its digits in a balanced 4-bit code, whose first bit splits 0-4 from
5-9 as a binary number's bits halve its values. IP addresses, versions and dotted dates stay integers; anything
unusual ('+', exponents, thousands separators, more than 18 digits) stays text around numbers, so that every way
of writing comes back exactly. Knowing the layout first, the model knows each digit's weight.

**Predicting each number from its role.** Slots remember the last values of a number by its place: the k-th
number of its line after the same text, the same JSON key, the same column of separators, the same local template.
Predictions: last value, trend, step, the line above, the number before it in the line and the same offset from
it, adaptive linear predictors, a median, a number of the line it copies (the title "M 1.2" copies the magnitude,
learned per slot). Each digit is coded with where each prediction stands against the digits still possible,
and each bit of its code against the value it splits, in eighths of the digit's weight.

**Computed numbers.** Programs print binary floats with more digits than they hold: Apple's prices are float32
printed with 6 decimals (`127.489998` is 127.49f), USGS coordinates float32 printed with 14 or 15 digits;
predicting the print of the nearest float32 saves 2% to 4%. A count of seconds since 1970 is often written again
as a date and a time in the same line: computed (the time zone learned by a vote that follows the change to winter
time), they become copies (bgl.log -12%). Bits that the bytes possible at a place decide are not coded.

| File | xz -9e | zpaq -m5 | paq8px -8 | gmz before | gmz |
|---|---|---|---|---|---|
| aapl.csv (stock prices) | 17,776 | 14,708 | 10,989 | 13,146 | **10,627** |
| intel.txt (sensor readings) | 164,824 | 116,944 | 74,290 | 91,849 | **73,072** |
| co2.txt (Mauna Loa) | 28,356 | 20,219 | 12,908 | 13,800 | 12,984 |
| meteo.csv (hourly weather) | 167,308 | 117,631 | 76,327 | 91,399 | 77,656 |
| bgl.log | 39,324 | 27,866 | 21,504 | 25,945 | 21,910 |
| hdfs.log | 42,348 | 35,425 | 25,783 | 30,447 | 26,896 |
| quakes.csv (USGS) | 235,796 | 176,454 | 134,938 | 174,751 | 144,678 |
| quakes.json (USGS) | 83,432 | 59,289 | 36,774 | 55,049 | 42,895 |
| openssh.log | 9,740 | 6,588 | 3,100 | 4,744 | 4,107 |
| apache.log | 6,744 | 4,742 | 2,488 | 3,145 | 3,046 |
| nci (33 MB, chemical structures) | | 1,251,218 | | 1,166,039 | **1,027,935** |

Over the ten files gmz goes from 504,275 to 417,871 bytes (-17.1%). It beats paq8px on two and comes within 2% on
three, about 10 times faster; prose does not change. What paq8px still wins on is the text around the numbers:
with every number replaced by #, quakes.json costs gmz 13,347 bytes against 7,075 for paq8px, place names
taking 41% of them.

Tried and dropped:
- residuals as sign and magnitude against a column prediction: worse everywhere but co2 (meteo +6.8%,
  quakes.json +15%);
- the value in binary rather than decimal digits: decimal wins on aapl (-17%), quakes and most logs, binary only
  on meteo;
- the text history without the digits: +0.1% to +0.9%;
- long copies from old lines (co2's "a year ago" column is not the value 52 lines above).

Several predictions are near the floor: the microseconds of bgl, the block identifiers of hdfs and the times of
earthquakes are essentially random.

The cost is speed: numeric text compresses 1.5 to 2.5 times slower than gmz's text model, and text just above 3%
of digits pays the comparison of two more plans for little (bib -1%, trans 0%).

Integrating it found that decoding a corrupted message could crash: a digit's first bits, coded as near certain,
came out otherwise, and a shift went negative. 126 decodings of corrupted messages under AddressSanitizer now all
end in a clean refusal.

Validated on 30 more files from new sources (finance, astronomy, genomics, logs, GPS, OpenStreetMap, government
CSV, JSON APIs, SQL dumps, C tables, and negative controls), the numeric plan is chosen on 23 and takes gmz 12.1%
below its text model overall, never worse, the controls untouched:

| Kind of data | Gain | Against paq8px |
|---|---|---|
| finance (Binance, EUR/USD ticks, ECB, VIX, CoinGecko) | 9% to 32% | -1.8% to +1.4% (ETHBTC trades +6.4%) |
| scientific tables (exoplanets, UCSC genes) | 24% to 29% | +3% to +5% |
| census CSV | 18% | -3.3% |
| SQL dumps, C tables of floats | 8% to 30% | +7% to +15% |
| genomics (VCF, GFF), geography (OSM, GPX, GeoJSON) | 5% to 14% | +4% to +27%, the text around the numbers |
| web logs (nginx, Zeek) | 0% to 0.6% | the plan rightly rejected |

Its failures: numbers that are categories (addresses, status codes, ports), random identifiers and hexadecimal
cut into decimal fragments, JSON on a single line (no line above to predict from), smooth functions (tables of
cosines), and columns computed from others. That last one is the axis again: in ETHBTC trades quote = price ×
quantity, in the census NPOPCHG = POP(t) - POP(t-1); a prototype learning, per column, a sum, difference or
product of two earlier numbers of the line takes ETHBTC 34% and the census 22.5% below paq8px.

### A million random digits

The file of Mark Nelson's challenge (AMillionRandomDigits.bin, 415,241 bytes, unbeaten since 2002) is the digits
of RAND's 1955 table read as one decimal number and written in binary: parsing the table's million digits into a
gmprime gives the same bytes, highest first. 10^6 log2(10) / 8 = 415,241.01 bytes, so the file holds exactly the
information of a million decimal digits, and any saving has to come from a flaw of the generator. Measured on the
digits, as 4-bit symbols in a gmprime:

| Test | Result | Random digits give |
|---|---|---|
| digit counts | chi-square 13.3 for 9 degrees of freedom | 9 (21.7 at 99%) |
| conditional entropy, orders 1 to 3 (Miller-Madow corrected) | 3.32192, 3.32191, 3.32179 bits | 3.32193 |
| adaptive code length, orders 1 to 3 | +127, +963, +6,523 bytes over the file | the cost of learning nothing |
| position in the 5-digit group, in the 50-digit line | chi-square 41.6 for 45, 454 for 450 | the degrees of freedom |
| longest repeated string | 11 digits | about 12 |

The exact enumerative bound: with the counts known, the digits are one of 2^3,321,834 arrangements, 95 bits fewer
than 10^6 digits take, but the counts cost 179 bits: +84 bits overall. `multiset_rank` computes that rank of a
million digits in 105 s. gmz gives 416,531 bytes (+1,290) and xz 415,328 (+87). No structure is visible with
these tools: a calibration, they find none where there is none.

## Speed

Starting point: about 150 KB/s per pass, and gmz compresses twice (with and without capitals, or full and bitmap
variants). Profile on lcet10.txt: `log` (stretch) 27%, table lookups and mixing 32%, learning 25%, `exp` (squash)
6%, copying mixer inputs 4%; at startup, fold analysis of the whole file and building the tables.

Times are for compressing in one pass (`-C off`) lcet10.txt, and with gmz's own choices ptt5 and kennedy.xls, on
Apple Silicon. Sizes change only where noted.

| Step | lcet10.txt | ptt5 | kennedy.xls | Size | Decision |
|---|---|---|---|---|---|
| start | 2.9 s | 5.3 s | 10.4 s | | |
| stretch and squash from tables (12 bits) | 2.4 s | 4.6 s | 9.1 s | +0.03% text, +1.4% and +3.1% on near-certain data | refined |
| ... stretch on 16 bits, squash over [-12, 12] | 2.3 s | 4.7 s | 9.0 s | within 2 bytes | kept |
| mixers read shared inputs; counter rates from a table | 2.3 s | 4.5 s | 8.6 s | identical | kept |
| tables from calloc; probability maps filled by copy | 2.3 s | 4.4 s | 8.6 s | identical | kept (less memory touched) |
| dot product in 4 partial sums | 2.2 s | 4.5 s | 8.5 s | identical | kept |
| slots grouped by nibble in 16-slot buckets | 1.8 s | 4.4 s | 8.1 s | +0.35% text, +2% kennedy.xls | candidate for a fast mode |
| prefetching the next bit's slots instead | 1.9 s | 4.4 s | 8.1 s | identical | kept |
| choosing variants on the first 64 KB, then one pass | 1.8 s* | 1.6 s | 5.9 s | +0.1% on 2 of 12 files | kept as the default; exhaustive on request |
| mixer weights and inputs in float instead of double | 1.3 s | 1.7 s | 5.6 s | same sizes (the messages differ) | kept |
| bit history updates from a transition table | 1.4 s | 1.7 s | 5.7 s | identical | dropped: no gain, the cost is the memory written, not the branches |

\* lcet10.txt with gmz's own choice of capitals: 3.1 s before, two passes.

Fold analysis on the first 128 KB instead of the whole file finds the same folds on all 12 files, 0.4 to 1 s
sooner.

Folded data, after the integer arithmetic, with identical messages: bit templates compiled once into runs of
consecutive distances, each read as one word (`compiled_template`, `stream_bits`), the aligned byte templates
prepared for each bit position, and pixel runs counted 56 bits at a time. kennedy.xls compresses 25% faster and
decompresses 27% faster, ptt5 16% and 26% (measured on a loaded machine, so as ratios). Prefetching the fold
contexts, or the first bit of each byte, gained nothing measurable: dropped.

After these steps the time is spread over the context lookups (33%), learning (25%), counter updates (29%) and
the mixers (9%): memory accesses rather than arithmetic.

### Result

gmz with everything kept, compared with the first version (compress / decompress):

| File | Before | After | Speed-up |
|---|---|---|---|
| alice29.txt | 2.75 s / 1.31 s | 1.16 s / 0.51 s | 2.4x / 2.6x |
| lcet10.txt | 7.28 s / 3.46 s | 1.84 s / 1.33 s | 4.0x / 2.6x |
| ptt5 | 7.81 s / 1.65 s | 1.79 s / 1.15 s | 4.4x / 1.4x |
| kennedy.xls | 15.6 s / 11.7 s | 6.0 s / 5.2 s | 2.6x / 2.3x |

Sizes moved by at most 2 bytes, except lcet10.txt (+95 bytes) and plrabn12.txt (+77 bytes), where the 64 KB
sample picks the capitals transform that the whole file does not.

After the integer arithmetic and the word reads of the fold templates (processor time, the machine being
loaded by other runs): alice29.txt 1.09 s / 0.48 s, lcet10.txt 1.86 s / 1.24 s, ptt5 1.26 s / 0.78 s,
kennedy.xls 4.62 s / 4.07 s, so 2.5 to 6.2 times faster than the first version to compress.

### Faster levels

Each lever measured alone over 8 files (floating point version, before the integer arithmetic), size change:

| Cost | Levers |
|---|---|
| under 0.1% | no line contexts (byte above), no last 3 bytes, no line-end distance, no class mixer, no capitals comparison (which compresses the sample twice) |
| 0.2% to 0.8% | nibble buckets (+0.32%), tables of 2^21 slots (+0.22%), orders 0 to 4 and 6 (+0.53%), no chained maps (+0.77%) |
| 1.4% to 3% | tables of 2^18 slots, orders 0 to 4, no word contexts, a single mixer |
| never | no prefetching (60% slower), no learning-rate schedule (larger and slower) |

Combined, the curve bends after dropping the chained maps: further levers cost 1.7% for 5% of time, then 2.8%
for 11%. Two levels were kept, stored in the message header (2 bits):

- `-1`: orders 0 to 4 and 6 in nibble buckets, without the line contexts, the last 3 bytes and the class
  mixer, tables of 2^21 slots, no capitals comparison;
- `-2`: the same without the chained maps.

Applied as is to folded data, `-1` cost kennedy.xls 13.7% for 1.4 times the speed: its byte contexts need orders
5 and 8 (+10% without them). Folded data now only leaves out its text contexts (words, lines, last
occurrences), +1.7% for the same speed.

| Level | Canterbury total | alice29.txt | lcet10.txt | plrabn12.txt | kennedy.xls |
|---|---|---|---|---|---|
| default | 333,834 | 1.09 s / 0.48 s | 1.86 s / 1.24 s | 1.99 s / 1.39 s | 4.62 s / 4.07 s |
| `-1` | 337,489 (+1.09%) | 0.43 s / 0.29 s | 0.91 s / 0.79 s | 1.00 s / 0.88 s | 3.22 s / 2.92 s |
| `-2` | 339,839 (+1.80%) | 0.39 s / 0.24 s | 0.87 s / 0.67 s | 0.88 s / 0.75 s | 2.98 s / 2.62 s |

## Portability: integer arithmetic

Floating point broke the decoding across platforms. The same gmz source, compressing the same files, gave 2 to 4
different messages per file over macOS and Linux, GCC and Clang, x86_64 and ARM64. Compilers may contract
`a * b + c` into a fused multiply-add where the processor has one, rounding once instead of twice, and do so
differently; `std::log` and `std::exp` differ between libraries too. One prediction off by the last bit and the
decoder goes astray for the rest of the message, silently. Sizes and times measured on one platform stay valid:
each build is consistent with itself.

The toolkit now computes on integers only, as lpaq does:

- probabilities on 16 bits (counters keep 22, so that the slow rates of frequent contexts still move them);
- the logistic domain in 1/256 units over [-12, 12]: `squash` interpolates 193 constant points, `stretch` is a
  table built from it;
- mixer weights in 8.24 fixed point, rates scaled by 65536, dot products on 64 bits;
- right shifts of negative numbers checked to be arithmetic at compile time.

Messages are now the same bytes on every build, and each decodes the others: on 6 Canterbury files, macOS (Apple
Clang) and Linux (GCC and Clang, x86_64 and ARM64) give one message per file. gmz also ends each message with a
CRC-32 of the original (4 bytes), so that a decoding gone astray fails instead of writing wrong data.

Without the checksum, sizes moved by at most 5 bytes (kennedy.xls 12,970 to 12,965, alice29.txt 35,563 to
35,562), and speed is the same or up to 16% better.

## Open ideas

Direction, after the broad evaluation: gmz's text model is the PAQ family's context mixing, the state of the art,
and adding to it would rebuild paq8px for a few percent. The results that broke ground came from recognizing the
computation behind the bytes: decimals stored as doubles (sao, 47% below zpaq), pixels as numbers (x-ray, below
paq8px), geometry found by folding the whole number. That is the axis now:

- Numbers written in text: done above. The structure text around them (names, messages) is where paq8px leads.
- Computed fields, which cost nothing once recognized: checksums (CRC of PNG chunks or zip entries), derived
  columns, counters, timestamps.
- Other computations to invert: floating point values from other operations than k * step, data already
  compressed (deflate inside PNG, docx, PDF).

Earlier ideas:

- Several folds at once, for data with nested periods.
- Records of variable length (osdb, where zpaq stays 0.4% ahead): resynchronizing on separators.
- An exact formula for sao's residuals: about 2.5 bits per double are left on the table.
- Elastic folds: rows of varying length (text lines, variable records), aligned on a separator.
- Linear complexity (Berlekamp-Massey) of bit planes and columns: the bits a linear recurrence would generate.
- Grammars (Re-Pair): the data as the output of a small program, the other way to treat a number as the result
  of a computation.
