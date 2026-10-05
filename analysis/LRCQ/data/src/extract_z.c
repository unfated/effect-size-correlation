/* Fast Pan-UKB flat-file parser for LRCQ (replaces the zcat | cut | awk pipe).
 * Usage: curl ... | extract_z <snp_universe.txt> | gzip > trait.z.gz
 * Reads a bgzipped Pan-UKB TSV on stdin, finds columns by header name
 * (af_EUR or af_controls_EUR, beta_EUR, se_EUR, low_confidence_EUR), and
 * prints one line per universe SNP ("chr:pos:ref:alt" order of the universe
 * file): the EUR Z-score, or NA if missing / low confidence / AF outside
 * [0.01, 0.99]. Exit code 2 if the header lacks a needed column.
 * Build: gcc -O2 -o extract_z extract_z.c -lz
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct { char *key; int idx; } Slot;
static Slot *tab; static size_t cap;
static unsigned long hsh(const char *s, size_t n) {
  unsigned long h = 1469598103934665603UL;
  for (size_t i = 0; i < n; i++) { h ^= (unsigned char)s[i]; h *= 1099511628211UL; }
  return h;
}
static void put(char *k, int idx) {
  size_t n = strlen(k), i = hsh(k, n) & (cap - 1);
  while (tab[i].key) i = (i + 1) & (cap - 1);
  tab[i].key = k; tab[i].idx = idx;
}
static int get(const char *k, size_t n) {
  size_t i = hsh(k, n) & (cap - 1);
  while (tab[i].key) {
    if (strlen(tab[i].key) == n && memcmp(tab[i].key, k, n) == 0) return tab[i].idx;
    i = (i + 1) & (cap - 1);
  }
  return -1;
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: extract_z universe.txt < sumstats.bgz\n"); return 1; }
  FILE *u = fopen(argv[1], "r"); if (!u) { perror("universe"); return 1; }
  size_t m = 0, ucap = 1 << 20; char **keys = malloc(ucap * sizeof(char *)); char line[1 << 16];
  while (fgets(line, sizeof line, u)) {
    line[strcspn(line, "\r\n")] = 0;
    if (m == ucap) { ucap *= 2; keys = realloc(keys, ucap * sizeof(char *)); }
    keys[m++] = strdup(line);
  }
  fclose(u);
  cap = 1; while (cap < 3 * m) cap <<= 1;
  tab = calloc(cap, sizeof(Slot));
  for (size_t i = 0; i < m; i++) put(keys[i], (int)i);
  double *z = malloc(m * sizeof(double)); char *have = calloc(m, 1);

  gzFile g = gzdopen(0, "rb"); gzbuffer(g, 1 << 20);
  char *buf = malloc(1 << 20);
  if (!gzgets(g, buf, 1 << 20)) { fprintf(stderr, "empty input\n"); return 1; }
  buf[strcspn(buf, "\r\n")] = 0;
  int c_af = -1, c_afc = -1, c_b = -1, c_s = -1, c_q = -1, col = 0;
  for (char *t = strtok(buf, "\t"); t; t = strtok(NULL, "\t"), col++) {
    if (!strcmp(t, "af_EUR")) c_af = col;
    else if (!strcmp(t, "af_controls_EUR")) c_afc = col;
    else if (!strcmp(t, "beta_EUR")) c_b = col;
    else if (!strcmp(t, "se_EUR")) c_s = col;
    else if (!strcmp(t, "low_confidence_EUR")) c_q = col;
  }
  if (c_af < 0) c_af = c_afc;
  if (c_af < 0 || c_b < 0 || c_s < 0 || c_q < 0) { fprintf(stderr, "missing EUR columns\n"); return 2; }
  int maxc = c_af; if (c_b > maxc) maxc = c_b; if (c_s > maxc) maxc = c_s; if (c_q > maxc) maxc = c_q;
  char *f[256]; char key[512]; long rows = 0;
  while (gzgets(g, buf, 1 << 20)) {
    rows++;
    int n = 0; char *p = buf; f[n++] = p;
    while (*p && *p != '\n' && n <= maxc) { if (*p == '\t') { *p = 0; f[n++] = p + 1; } p++; }
    if (n <= maxc) continue;
    char *e = strchr(f[maxc], '\t'); if (e) *e = 0; e = strchr(f[maxc], '\n'); if (e) *e = 0;
    int kl = snprintf(key, sizeof key, "%s:%s:%s:%s", f[0], f[1], f[2], f[3]);
    int idx = get(key, (size_t)kl); if (idx < 0) continue;
    if (strcmp(f[c_q], "false") || !strcmp(f[c_b], "NA") || !strcmp(f[c_s], "NA") || !strcmp(f[c_af], "NA")) continue;
    double af = atof(f[c_af]), b = atof(f[c_b]), s = atof(f[c_s]);
    if (s <= 0 || af < 0.01 || af > 0.99) continue;
    z[idx] = b / s; have[idx] = 1;
  }
  int err; const char *msg = gzerror(g, &err);
  if (err != Z_OK && err != Z_STREAM_END) { fprintf(stderr, "gz error: %s\n", msg); return 3; }
  for (size_t i = 0; i < m; i++) { if (have[i]) printf("%.5g\n", z[i]); else fputs("NA\n", stdout); }
  fprintf(stderr, "rows %ld\n", rows);
  return 0;
}
