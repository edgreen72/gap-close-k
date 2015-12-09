#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define MAX_FN_LEN (1024)
#define MAX_SEQ_LEN (8388608) // ~8 Mb
#define MAX_ID_LEN (256)
#define MAX_KMER_OCCUR (128)
#define MAX_UNIQ_KMER (3)
#define K (64)
#define NUM_N (5)
#define FLANK_SIZE (100)
#define MIN_SPANNERS (2)
#define KMER_DEF (14)
#define MAX_NUM_SPANNERS (2048) // 2k spanners
#define MAX_FQ_LEN (1024) // maximum length of a sequenc read in fastq file
#define CONS_LENGTH_CUT (95) // minimum percent of spanning reads to consensus to close gap
#define CONS_BASE_PERC (90) // minimum percent of bases that must match to call a consensus base at each position
#define DEBUG (0)
#define VERSION (6)

int debug_info = DEBUG; // a true global variable - Use this in functions directly, i.e. not passed


typedef struct kmers {
  unsigned int k; // length of the kmers
  short unsigned int* ka; // the kmer array
} Kmers;
typedef struct kmers* KSP;

typedef struct gaps {
  size_t* gap_starts;
  size_t* gap_ends;
  size_t num_gaps;
} Gaps;
typedef struct gaps* GP;
typedef struct spanners {
  size_t num_spanners;
  char* spanners[MAX_NUM_SPANNERS];
  char span_seq[ MAX_NUM_SPANNERS * (MAX_FQ_LEN + 1) ];
  char cons_seq[ MAX_FQ_LEN + 1];
  int As[MAX_FQ_LEN+1];
  int Cs[MAX_FQ_LEN+1];
  int Gs[MAX_FQ_LEN+1];
  int Ts[MAX_FQ_LEN+1];
} Spanners;
typedef struct spanners* SP;
typedef struct fastq_read {
  char seq[MAX_FQ_LEN+1];
  char qual[MAX_FQ_LEN+1];
  char seqrc[MAX_FQ_LEN+1];
  char qualrc[MAX_FQ_LEN+1];
  char id[MAX_ID_LEN+1];
} Fastq_read;
typedef struct fastq_read* FQRP;
void make_gap_table(  const char* genome_fn, const KSP ks, const int num_n,
                      const int flank_size );
void fasta_close_gaps( const char* genome_fn, const char* fastq_fn,
                       const KSP ks, const int num_n, const int flank_size,
                       const int min_spanners );
void populate_spans( SP spans, const char k1[], const char k2[], const char fastq_fn[] );
size_t consensus_length( const SP spans );
void find_cons_seq( SP spans, size_t cons_length );
char best_base_from_counts( int* perc, const int As, const int Cs, const int Gs, const int Ts );
int close_unique_k( size_t beg, size_t end, const char* seq, const KSP ks, char* kmer,
		    size_t* mko );
void populate_GP( const char* seq, GP gaps, const int num_n );
SP init_Spans( void );
void reset_spans( SP spans );
GP init_Gaps( size_t seq_len, const int num_n );
void free_Gaps( GP gaps );
int fasta_pop_kmers( const char genome_fn[], KSP ks );
size_t next_fa( FILE* fp, char* seq, char* id );
KSP init_KSP( int k );
FILE * fileOpen(const char *name, char access_mode[]);
int kmer2inx( const char* kmer,
              const unsigned int kmer_len,
              size_t* inx );
void output_kmer_dist( KSP ks );
int next_fastq( FILE* fastq, FQRP fqrp );
inline char revcom_char(const char base);
void revcom_seq( const char seq[], char rcom_seq[] );
void rev_qual( const char q[], char qrc[] );

void help( void ) {
  printf( "read-close-gaps VERSION %d\n", VERSION );
  printf( "   -g <genome fasta file>\n" );
  printf( "   -f <fastq file>\n" );
  printf( "   -k <kmer k; default = %d\n", KMER_DEF );
  printf( "   -n <number of consecutive Ns to define a gap; default = %d\n", NUM_N );
  printf( "   -F <size of Flank region to find unique kmers; default = %d\n", FLANK_SIZE );
  printf( "   -m <minimum kmer-containing, gap-spanning reads to close a gap; default = %d\n", MIN_SPANNERS );
  printf( "   -t <send a table to STDOUT with ID GSTART GEND K1POS K2POS K1 K2>\n" );
  printf( "   -d <DEBUG mode - print a bunch of info to STDERR along the way>\n" );
  printf( "If the -t option is given, read-close-gaps will not attempt to close any\n" );
  printf( "gaps. Instead, it will just make a table of the gaps seen in the input\n" );
  printf( "genome and the kmers that it would attempt to use to close these gaps.\n" );
  exit( 0 );
  
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  int k = KMER_DEF;
  int flank_size = FLANK_SIZE;
  int num_n = NUM_N;
  int min_spanners = MIN_SPANNERS;
  int ich;
  char genome_fn[MAX_FN_LEN];
  char fastq_fn[MAX_FN_LEN];
  KSP ks;
  int make_GT = 0;

  /*          PROCESS ARGUMENTS            */
  if ( argc == 1 ) {
    help();
  }
  
  while( (ich=getopt( argc, argv, "g:k:n:f:F:m:td" )) != -1 ) {
    switch(ich) {
    case 'g' :
      strcpy( genome_fn, optarg );
      break;
    case 'k' :
      k = atoi( optarg );
      break;
    case 'n' : // Number of consective Ns to define a gap
      num_n = atoi( optarg );
      break;
    case 'F' :
      flank_size = atoi( optarg ); // Size region surrounding gaps to look for unique kmers
      break;
    case 'f' :
      strcpy( fastq_fn, optarg );
      break;
    case 'm' : // Minimum number of reads spanning gap to close it
      min_spanners = atoi( optarg );
      break;
    case 't' :
      make_GT = 1;
      break;
    case 'd' :
      debug_info = 1;
      break;
    default :
      help();
    }
  }

  /*          INITIALIZE            */
  ks = init_KSP( k );
  fprintf(stderr,  "Initialized kmer structure\n" );
  /* Go through the input genome, one fasta record at a time */
  fasta_pop_kmers( genome_fn, ks );
  fprintf(stderr,  "Scanned genome for kmers k = %d\n", ks->k );
  
  //output_kmer_dist( ks );

  if ( make_GT ) {
    make_gap_table( genome_fn, ks, num_n, flank_size );
  }
  else {
    /* Go through the input genome, again, one fasta record at a time */
    fasta_close_gaps( genome_fn, fastq_fn, ks, num_n,
                      flank_size, min_spanners );

  }
  exit( 0 );
}

/* make_gap_table */
void make_gap_table(  const char* genome_fn, const KSP ks, const int num_n,
                      const int flank_size ) {
  FILE* fa;
  GP gaps_p;
  SP spanners;
  size_t seq_len, beg, end, i, k1_mko, k2_mko;
  int k1_pos, k2_pos;
  char k1[K+1]; // kmer upstream of gap
  char k2[K+1]; // kmer downstream of gap
  char id[MAX_ID_LEN+1];
  char* seq;

  if ( debug_info ) {
    fprintf( stderr, "Starting to make gap table\n" );
    fflush(stderr);
  }

  fa = fileOpen( genome_fn, "r" );
  seq = (char*)malloc(sizeof(char)*(MAX_SEQ_LEN+1));

  while( seq_len = next_fa( fa, seq, id ) ) {
    if ( debug_info ) {
      fprintf( stderr, "Scanning %s for gaps\n", id );
      fflush(stderr);
    }

    /* Initialize the structure for keeping gap positions */
    gaps_p = init_Gaps( seq_len, num_n );

    /* Scan the sequence for runs of Ns and populate the gaps_p */
    populate_GP( seq, gaps_p, num_n );

    /* Find the k1 and k2 around each gap */
    for( i = 0; i < gaps_p->num_gaps; i++ ) {
      if ( gaps_p->gap_starts[i] < flank_size ) {
        beg = 0;
      }
      else {
        beg = gaps_p->gap_starts[i] - flank_size;
      }
      k1_pos = close_unique_k( beg, gaps_p->gap_starts[i], seq, ks, k1, &k1_mko );

      if ( gaps_p->gap_ends[i] + flank_size > seq_len ) {
        end = seq_len;
      }
      else {
        end = gaps_p->gap_ends[i] + FLANK_SIZE;
      }
      k2_pos = close_unique_k( end, gaps_p->gap_ends[i], seq, ks, k2, &k2_mko );

      /* If we found suitable k1 and k2, then just print some info on them */
      if ( (k1_pos >= 0) &&
           (k2_pos >= 0) ) {
        printf( "%s %d %d %d %d %s %s %d %d\n",
                id, (int)gaps_p->gap_starts[i], (int)gaps_p->gap_ends[i],
                (int)k1_pos, (int)k2_pos, k1, k2, (int)k1_mko, (int)k2_mko );
      }
      else {
        printf( "%s %d %d NA NA NA NA NA NA\n",
                id, (int)gaps_p->gap_starts[i], (int)gaps_p->gap_ends[i] );
      }
    }
    free_Gaps( gaps_p ); // free the memory we had allocated for this GP and its guts
  }
}


/* fasta_close_gaps
   ARGS: (3) const char genome_fn - the filename of the fasta format genome sequence
             const char fastq_fn  - the filename of the fastq format reads that might
                                    span gaps
             const KSP ks         - the kmers structure 
   RETURNS: nothing (for now)
   Goes through each sequence in the genome and finds gaps (runs of Ns). Then, finds
   the most unique and closest kmers flanking the gap. Then, finds reads containing
   both of those kmers and determines if they are unanimous about the distance between
   the kmers. If so, take the consensus sequence from the reads that span the gap to
   close the gap.
*/
void fasta_close_gaps( const char* genome_fn, const char* fastq_fn,
                       const KSP ks, const int num_n, const int flank_size,
                       const int min_spanners ) {
  FILE* fa;
  FILE* fq;
  GP gaps_p;
  SP spanners;
  size_t seq_len, beg, end, cons_length, i, k1_mko, k2_mko;
  int k1_pos, k2_pos;
  SP spans;
  char k1[K+1]; // kmer upstream of gap
  char k2[K+1]; // kmer downstream of gap
  char id[MAX_ID_LEN+1];
  char* gen_gap_seq;
  char* seq;

  if ( debug_info ) {
    fprintf( stderr, "Starting gap closing\n" );
    fflush(stderr);
  }

  fa = fileOpen( genome_fn, "r" );
  fq = fileOpen( fastq_fn, "r" );
  spans = init_Spans();
  spans->num_spanners = 0; // unnecessary
  seq = (char*)malloc(sizeof(char)*(MAX_SEQ_LEN+1));

  while( seq_len = next_fa( fa, seq, id ) ) {
    if ( debug_info ) {
      fprintf( stderr, "Scanning %s for gaps\n", id );
      fflush(stderr);
    }

    /* Initialize the structure for keeping gap positions */
    gaps_p = init_Gaps( seq_len, num_n );

    /* Scan the sequence for runs of Ns and populate the gaps_p */
    populate_GP( seq, gaps_p, num_n );

    /* Try to close each gap */
    for( i = 0; i < gaps_p->num_gaps; i++ ) {
      if ( gaps_p->gap_starts[i] < flank_size ) {
        beg = 0;
      }
      else {
        beg = gaps_p->gap_starts[i] - flank_size;
      }
      k1_pos = close_unique_k( beg, gaps_p->gap_starts[i], seq, ks, k1, &k1_mko );

      if ( gaps_p->gap_ends[i] + flank_size > seq_len ) {
        end = seq_len;
      }
      else {
        end = gaps_p->gap_ends[i] + FLANK_SIZE;
      }
      k2_pos = close_unique_k( end, gaps_p->gap_ends[i], seq, ks, k2, &k2_mko );

      /* If we found suitable k1 and k2, then search the reads for them */
      if ( (k1_pos >= 0) &&
           (k2_pos >= 0) ) {
        reset_spans( spans );
        populate_spans( spans, k1, k2, fastq_fn );
        if ( spans->num_spanners < min_spanners ) {
          printf( "Not close gap in %s between %d and %d - only %d spanning reads\n",
                  id, k1_pos, (k2_pos + ks->k), (int)spans->num_spanners );
        }
        else {
          cons_length = consensus_length(spans);
          if ( cons_length ) {
            gen_gap_seq = (char*)malloc(sizeof(char) * (k2_pos - k1_pos + ks->k + 1));
            strncpy( gen_gap_seq, &seq[k1_pos], (k2_pos - k1_pos + ks->k) );
            gen_gap_seq[k2_pos - k1_pos + ks->k] = '\0';
            find_cons_seq( spans, cons_length );
            printf( "Close gap in %s between %d and %d with %d spanning reads. %d and %d kmer occurances\n", 
                    id, (int)k1_pos, (int)(k2_pos + ks->k), (int)spans->num_spanners,
		    (int)k1_mko, (int)k2_mko );
            printf( "Genome: %s\nConsen: %s\n\n",
                    gen_gap_seq,
                    spans->cons_seq );
            free( gen_gap_seq );
          }
          else {
            printf( "Did not close gap in %s between %d and %d. No consensus length\n",
                    id, (int)gaps_p->gap_starts[i], (int)gaps_p->gap_ends[i] ); 
          }
        }
      }
      else { // didn't find one or more of k1, k2
        printf( "Could not attempt a gap closure in %s between %d and %d\n",
                id, (int)gaps_p->gap_starts[i], (int)gaps_p->gap_ends[i] );
      }
    }
    free_Gaps( gaps_p ); // free the memory we had allocated for this GP and its guts
  }
}

/* reset_spans */
void reset_spans( SP spans ) {
  size_t i;
  spans->num_spanners = 0;
  for( i = 0; i < (MAX_FQ_LEN + 1); i++ ) { // could use memset
    spans->As[i] = 0;
    spans->Cs[i] = 0;
    spans->Gs[i] = 0;
    spans->Ts[i] = 0;
  }
  return;
}

/* populate_spans
   ARGS (4) SP spans - the spanners data structure that will hold the kmer spanning
                       sequences we find
            const char k1[] - the first kmer
            const char k2[] - the second kmer
            const char fastq_fn[] - the filename of the fastq reads to search
    RETURNS - nothing
    Searches the reads in the specified fastq file, forward and revcom, and populates
    spans with the sequence segments that contain k1 and k2 in the proper order
 */

void populate_spans( SP spans, const char k1[], const char k2[], const char fastq_fn[] ) {
  FILE* fq;
  FQRP fqrp; // pointer to structure to hold a fastq read (forward and reverse)
  char* foundk1;
  char* foundk2;
  char k1tok2[MAX_FQ_LEN+1];
  size_t len;

  fq = fileOpen( fastq_fn, "r" );
  fqrp = (FQRP)malloc(sizeof(struct fastq_read));

  while( next_fastq( fq, fqrp ) == 1 ) {
    revcom_seq( fqrp->seq, fqrp->seqrc );

    foundk1 = strstr( fqrp->seq, k1 );
    if ( foundk1 == NULL ) {
      ;
    }
    else {
      foundk2 = strstr( foundk1, k2 );
      if ( foundk2 == NULL ) {
        ;
      }
      else {
        len = (size_t)(foundk2-foundk1 + strlen(k2));
        strncpy( k1tok2, foundk1, len );
        k1tok2[len] = '\0';
        if ( spans->num_spanners < MAX_NUM_SPANNERS ) {
          strcpy( spans->spanners[spans->num_spanners], k1tok2 );
          spans->num_spanners++;
        }
        else { // we've filled up spanners
          fclose( fq );
          free(fqrp);
          return;
        }
      }
    }

    /* Now check revcom */
    foundk1 = strstr( fqrp->seqrc, k1 );
    if ( foundk1 == NULL ) {
      ;
    }
    else {
      foundk2 = strstr( foundk1, k2 );
      if ( foundk2 == NULL ) {
        ;
      }
      else {
        len = (size_t)(foundk2-foundk1 + strlen(k2));
        strncpy( k1tok2, foundk1, len );
        k1tok2[len] = '\0';
        if ( spans->num_spanners < MAX_NUM_SPANNERS ) {
          strcpy( spans->spanners[spans->num_spanners], k1tok2 );
          spans->num_spanners++;
        }
        else { // we've filled up spanners
          fclose( fq );
          free(fqrp);
          return;
        }
      }
    }
    
  }
  fclose( fq );
  free(fqrp);
  return;
}

size_t consensus_length( const SP spans ) {
  size_t i;
  size_t max_count = 0;
  size_t max_len;
  size_t span_lengths[ MAX_FQ_LEN + 1];
  float most_popular_perc;
  /* reused memory - zero it out! */
  for( i = 0; i <= MAX_FQ_LEN; i++ ) {
    span_lengths[i] = 0;
  }
  /* Populate span_lengths with the counts for all lengths */
  for ( i = 0; i < spans->num_spanners; i++ ) {
    span_lengths[ strlen(spans->spanners[i]) ]++;
  }

  /* Find the most common length */
  for( i = 1; i < MAX_FQ_LEN; i++ ) {
    if ( span_lengths[i] > max_count ) {
      max_len = i;
      max_count = span_lengths[i];
    }
  }
  if ( max_count > 0 ) {
    most_popular_perc = (float)max_count/spans->num_spanners;
    most_popular_perc *= 100;
    if ( most_popular_perc >= CONS_LENGTH_CUT ) {
      return max_len;
    }
  }
  return 0;
}

void find_cons_seq( SP spans, size_t cons_length ) {
  size_t i, pos; 
  int perc;
  char best_base;
  for( i = 0; i < spans->num_spanners; i++ ) {
    if ( strlen( spans->spanners[i] ) == cons_length ) {
      for( pos = 0; pos < cons_length; pos++ ) {
        switch (spans->spanners[i][pos]) {
        case 'A' :
          spans->As[pos]++;
          break;
        case 'C' :
          spans->Cs[pos]++;
          break;
        case 'G' :
          spans->Gs[pos]++;
          break;
        case 'T' :
          spans->Ts[pos]++;
          break;
        }
      }
    }
    else { 
      ;
    }
  }
  for( pos = 0; pos < cons_length; pos++ ) {
    best_base = best_base_from_counts( &perc, 
                                       spans->As[pos], 
                                       spans->Cs[pos], 
                                       spans->Gs[pos], 
                                       spans->Ts[pos] );
    if ( perc >= CONS_BASE_PERC ) {
      spans->cons_seq[pos] = best_base;
    }
    else {
      spans->cons_seq[pos] = 'N';
    }
  }
  spans->cons_seq[pos] = '\0';
}

char best_base_from_counts( int* perc, const int As, const int Cs, const int Gs, const int Ts ) {
  int total = 0;
  total = As + Cs + Gs + Ts;
  if ( (As > Cs) &&
       (As > Gs) &&
       (As > Ts ) ) {
    *perc = (As / total) * 100;
    return 'A';
  }
  if ( (Cs > As) &&
       (Cs > Gs) &&
       (Cs > Ts) ) {
    *perc = (Cs / total) * 100;
    return 'C';
  }
  if ( (Gs > As) &&
       (Gs > Cs) &&
       (Gs > Ts) ) {
    *perc = (Gs / total) * 100;
    return 'G';
  }
  if ( (Ts > As) &&
       (Ts > Cs) &&
       (Ts > Gs) ) {
    *perc = (Ts / total) * 100;
    return 'T';
  }
  *perc = 25;
  return 'N';
}

/* close_unique_k
   ARGS (5) size_t beg - beginning coordinate of region to search
            size_t end - end coordinate of region to search
            const char* seq - the sequence to search
            const KSP ks - the kmer structure
            char* kmer - the found kmer
   RETURNS: position of most unique kmer, closest to gap; -1 if none found
   Searches the region of seq indicated by beg and end. Finds the most unique kmer
   in the region. In case of tie (common) will find the kmer that is both most
   unique and closest to the gap. If the first coordinate is smaller, the gap is
   downstream of the region mentioned. If the first coordinate is bigger, the gap
   is upstream of the region mentioned
*/
int close_unique_k( size_t beg, size_t end, const char* seq, const KSP ks,
		    char* kmer, size_t* mko ) {
  size_t min_kmer_occur = MAX_KMER_OCCUR; // initialize this; used to keep track of how many
  // times the kmer that is lowest in genome-wide occurance is seen
  size_t pos;
  size_t min_k_pos = 0; // keep the position of the most rare kmer found
  size_t inx;

  /* Check to make sure region is big enough for at least one k-mer */
  if ( abs( end-beg ) < ks->k ) {
    return -1;
  }
  
  if ( beg < end ) { // gap starts at end
    for( pos = beg; pos < (end - ks->k); pos++ ) {
      if ( kmer2inx( &seq[pos], ks->k, &inx ) ) {
        if ( ks->ka[inx] <= min_kmer_occur ) {
          min_kmer_occur = ks->ka[inx];
          min_k_pos = pos;
        }
      }
    }
  }
  else { // go backwards
    for( pos = (beg - ks->k); pos >= end; pos-- ) {
      if ( kmer2inx( &seq[pos], ks->k, &inx ) ) {
        if ( ks->ka[inx] <= min_kmer_occur ) {
          min_kmer_occur = ks->ka[inx];
          min_k_pos = pos;
        }
      }
    }
  }
  *mko = min_kmer_occur;
  if (min_kmer_occur < MAX_UNIQ_KMER ) { // find ANYTHING?
    strncpy( kmer, &seq[min_k_pos], ks->k );
    kmer[ks->k] = '\0';
    return (int)min_k_pos;
  }
  else {
    return -1;
  }
}

/* populate_GP
   ARGS: (1) const char* seq - the sequence of this chromosome/scaffold/whatever
             that might have some run of Ns in it that are gaps
         (2) GP gaps - the data structure to be populated
   RETURNS: nothing
   Scans the sequence in seq looking for runs of >= NUM_N Ns in a row.
   Each region found goes into the gaps structure with its beginning and end
   0-indexed, open-ended
*/
void populate_GP( const char* seq, GP gaps, const int num_n ) {
  size_t pos;
  size_t seq_len;
  size_t curr_start;
  seq_len = strlen( seq );
  int N_region = 0;
  for( pos = 0; pos <= seq_len; pos++ ) {
    if ( (seq[pos] == 'N') ||
         (seq[pos] == 'n') ) { // current position is an N
      if ( N_region ) {
        ; // just continuing a run of Ns
      }
      else {
        // starting a new run of Ns
        curr_start = pos;
      }
      N_region = 1;
    }
    else { // current position is not an N
      if ( N_region ) {
        /* just finished a run of Ns. Is it long enough
           to be considered a gap */
        if ( (pos - curr_start) > NUM_N ) {
          gaps->gap_starts[gaps->num_gaps] = curr_start;
          gaps->gap_ends[gaps->num_gaps]   = pos;
          gaps->num_gaps++;
        }
      }
      N_region = 0;
    }
  }
}

/* init_Spans
 */
SP init_Spans( void ) {
  SP spans;
  size_t i;
  spans = (SP)malloc(sizeof( struct spanners ));
  for( i = 0; i < MAX_NUM_SPANNERS; i++ ) {
    spans->spanners[i] = &(spans->span_seq[ i * (MAX_FQ_LEN + 1)]);
  }
  spans->num_spanners = 0;
  return spans;
}

/* init_Gaps
   ARGs: (1) size_t seq_len - the length of the sequence we'll scan for gaps
   positions
   RETURNS: GP (pointer to struct gaps
   This is dynamically allocated and freed by free_Gaps(GP)
*/
GP init_Gaps( size_t seq_len, const int num_n ) {
  GP gaps;
  gaps = (GP)malloc(sizeof( struct gaps ));
  /* The most gaps we COULD have is the length of this sequence divided by
     the number of Ns that define a gap, plus one */
  gaps->gap_starts = (size_t*)malloc( sizeof( size_t ) * (seq_len/(num_n+1)));
  gaps->gap_ends   = (size_t*)malloc( sizeof( size_t ) * (seq_len/(num_n+1)));
  gaps->num_gaps   = 0;
  return gaps;
}

/* free_Gaps
   ARGS: (1) pointer to a struct gaps
   RETURNS: nothing
   Frees all the memory allocated for this struct gaps, including
   internal memory
*/
void free_Gaps( GP gaps ) {
  free( gaps->gap_starts );
  free( gaps->gap_ends );
  free( gaps );
}

/* Return : 0 => copacetic
   non-zero => problem */
int fasta_pop_kmers( const char genome_fn[], KSP ks ) {
  char id[MAX_ID_LEN+1];
  char* seq;
  char c;
  int pos = 0;
  size_t inx = 0;
  size_t seq_len;
  FILE* fa;
  int k = ks->k;
  int kmer_status;

  fa = fileOpen( genome_fn, "r" );
  seq = (char*)malloc(sizeof(char)*(MAX_SEQ_LEN+1));

  /* Go through the fasta file, one sequence at a time */
  while( seq_len = next_fa( fa, seq, id ) ) {
    if ( debug_info ) {
      fprintf( stderr, "Counting kmers in %s\n", id );
    }
    if ( seq_len >= ks->k ) { // need at least one k-mer to count kmers!
      for( pos = 0; pos < seq_len - ks->k; pos++ ) {
        kmer_status = kmer2inx( &seq[pos], ks->k, &inx );
        if ( kmer_status ) {
          if ( ks->ka[inx] < MAX_KMER_OCCUR ) {
            ks->ka[inx]++;
          }
        }
      }
    }
  }
  free(seq);
  fclose(fa);
  return 0;
}

/* Takes a filehandle of a fasta file at the '>' character
   Reads the next sequence and puts it in the char* seq 
*/
size_t next_fa( FILE* fp, char* seq, char* id ) {
  char c;
  size_t i = 0;
  c = fgetc(fp);
  if ( c == '>' ) {
    c = fgetc(fp);
    while ( !isspace(c) ) { // load up the ID
      id[i++] = c;
      c = fgetc(fp);
    }
    id[i] = '\0';
    while( c != '\n' ) {
      c = fgetc(fp);
    }
    i = 0;
    while( (c != '>') &&
           (c != EOF) &&
           (i < MAX_SEQ_LEN) ) {
      if ( isspace(c) ) {
        ;
      }
      else {
        c = toupper(c);
        seq[i++] = c;
      }
      c = fgetc( fp );
    }
    seq[i] = '\0';
    ungetc(c, fp);
  }
  else {
    if ( c == EOF ) {
      return 0;
    }
  }
  if ( i == MAX_SEQ_LEN ) {
    fprintf( stderr, "%s is truncated to %d\n", id, MAX_SEQ_LEN );
  }
  return i;
}


/* Initialize a new KSP
   Populate all ka array elements with NULL
*/
KSP init_KSP( int k ) {
  KSP ks;
  size_t i, len;
  ks = (KSP)malloc(sizeof(Kmers));
  len = 1<<(k*2);
  ks->ka = (short unsigned int*)malloc(sizeof(short unsigned int)*len);

  ks->k = k;
  for( i = 0; i < len; i++ ) {
    ks->ka[i] = 0;
  }
  return ks;
}


/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}

/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
             might not be null-terminated
         (2) length of the kmer
         (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
              const unsigned int kmer_len,
              size_t* inx ) {
  size_t l_inx = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
    case 'A' :
      l_inx += 0;
      break;
    case 'C' :
      l_inx += 1;
      break;
    case 'G' :
      l_inx += 2;
      break;
    case 'T' :
      l_inx += 3;
      break;
    default :
      return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}

void output_kmer_dist( KSP ks ) {
  /* dist[] => index is the number of occurances for a kmer
               value is the number of kmers that occur that
               number of times */
  size_t dist[MAX_KMER_OCCUR + 1];
  size_t i;
  size_t num_kmers;
  num_kmers = 1<<(2 * ks->k);

  for( i = 0; i < (MAX_KMER_OCCUR + 1); i++ ) {
    dist[i] = 0;
  }
  for( i = 0; i < num_kmers; i++ ) {
    dist[ ks->ka[i]]++;
  }

  for( i = 0; i <= MAX_KMER_OCCUR; i++ ) {
    printf( "%d %d\n", (int)i, (int)dist[i] );
  }
  return;
}

/* next_fastq
 */
int next_fastq( FILE* fastq, FQRP fqrp ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    fqrp->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      fqrp->id[i] = '\0';
    }
  }
  fqrp->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      fqrp->seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  fqrp->seq[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", fqrp->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      fqrp->qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  fqrp->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return 0;
  }
  return 1;
}

inline char revcom_char(const char base) {
  switch (base) {
  case 'A':
    return 'T';
  case 'a' :
    return 't';

  case 'C':
    return 'G';
  case 'c' :
    return 'g';

  case 'G':
    return 'C';
  case 'g' :
    return 'c';
    
  case 'T':
    return 'A';
  case 't' :
    return 'a';

  case '-':
    return '-';
    
  case 'N':
    return 'N';
  case 'n':
    return 'n';

  case 'X':
    return 'X';
  case 'x':
    return 'x';
    
  default:
    fprintf( stderr, "Do not know how to revcom \"%c\"\n", base);
    return 'N';
  }
}

void revcom_seq( const char seq[], char rcom_seq[] ) {
  int len, i;
  len = strlen(seq);

  for (i = 0; i < len; i++) {
    rcom_seq[i] = revcom_char(seq[i]);
  }
}

void rev_qual( const char q[], char qrc[] ) {
  int len, i;
  len = strlen(q);
  
  for (i = 1; i <= len; i++) {
    qrc[len - i] = q[i-1];
  }
}
