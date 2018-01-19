#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define MAX_FN_LEN (1024)
#define MAX_SEQ_LEN (450000000) // 450 Mb
#define MAX_RECORDS (500000) // Maximum number of contigs/scaffolds/whatevers
#define MAX_ID_LEN (256)
#define DEBUG (0)
#define VERSION (1)

int base_comp[ 6 ]; // counts for each base: A, C, G, T, N, all others
int dinuc_comp[ 16 ]; // counts for each valid dinucleotide
int N_runs[10]; // counts for number of runs of Ns of various log10 lengths
size_t total_len = 0;
int total_records = 0;
size_t lengths[ MAX_RECORDS ]; // lengths of all fasta sequences
int debug_info = DEBUG; // a true global variable
// Use this in functions directly, i.e. not passed



int fasta_analyze( const char genome_fn[] );
size_t next_fa( FILE* fp, char* seq, char* id );
FILE * fileOpen(const char *name, char access_mode[]);
void report( void );
void increment_base_comp( const char b );
void increment_dinuc_comp( const char b, const char c );
static int cmpintp( const void *p1, const void *p2 );

void help( void ) {
  printf( "describe-assembly VERSION %d\n", VERSION );
  printf( "   -g <genome fasta file>\n" );
  printf( "   -d <DEBUG mode - print a bunch of info to STDERR along the way>\n" );
  printf( "Gives summary statistics for base composition, dinucleotide frequencies,\n" );
  printf( "N50, contig size distribution, number and sizes of runs of Ns.\n" );
  exit( 0 );
  
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  extern int optin;
  int ich;
  char genome_fn[MAX_FN_LEN];

  /*          PROCESS ARGUMENTS            */
  if ( argc == 1 ) {
    help();
  }
  
  while( (ich=getopt( argc, argv, "g:d" )) != -1 ) {
    switch(ich) {
    case 'g' :
      strcpy( genome_fn, optarg );
      break;
    default :
      help();
    }
  }

  /* Go through the input genome, one fasta record at a time */
  fasta_analyze( genome_fn );

  /* Make output */
  report();
  exit( 0 );
}

/* Returns size_t => length of fasta sequence
   0 => problem */
int fasta_analyze( const char genome_fn[] ) {
  char id[MAX_ID_LEN+1];
  char* seq;
  char c;
  int pos = 0;
  int N_region = 0;
  int N_len;
  size_t inx = 0;
  size_t seq_len;
  size_t i = 0;
  FILE* fa;

  fa = fileOpen( genome_fn, "r" );
  seq = (char*)malloc(sizeof(char)*(MAX_SEQ_LEN+1));

  /* Go through the fasta file, one sequence at a time */
  while( seq_len = next_fa( fa, seq, id ) ) {
    lengths[ i++ ] = seq_len;
    total_records++;
    total_len += seq_len;
    N_region = 0;
    if ( debug_info ) {
      fprintf( stderr, "Analyzing %s\n", id );
    }
    for( pos = 0; pos <= seq_len - 1; pos++ ) {
      if ( !isspace( seq[pos] ) ) {
	  increment_base_comp( seq[pos] );
	  if ( !isspace( seq[pos + 1] ) ) {
	    increment_dinuc_comp( seq[pos], seq[pos+1] );
	  }
      }
      
      if ( (seq[pos] == 'N') ||
	   (seq[pos] == 'n') ) {
	if ( N_region ) {
	  N_len++;
	}
	else {
	  N_len    = 1;
	  N_region = 1;
	}
      }
      else {
	if ( N_region ) {
	  N_runs[ (int)floor(log10( N_len )) ]++;
	  N_region = 0;
	}
      }
    }
  }
  free(seq);
  fclose(fa);
  return 0;
}


void increment_base_comp( const char b ) {
  switch(toupper(b)) {
  case 'A' :
    base_comp[0]++;
    break;
  case 'C' :
    base_comp[1]++;
    break;
  case 'G' :
    base_comp[2]++;
    break;
  case 'T' :
    base_comp[3]++;
    break;
  case 'N' :
    base_comp[4]++;
    break;
  default :
    base_comp[5]++;
  }
}

void increment_dinuc_comp( const char b, const char c ) {
  size_t l_inx = 0;
  switch(toupper(b)) {
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
    return; // not valid - just quit
  }
  l_inx = l_inx << 2;
  switch(toupper(c)) {
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
    return; // not valid - just quit
  }
  dinuc_comp[ l_inx ]++;
}

void report( void ) {
  size_t i;
  size_t n50_tmp = 0;
  printf( "# Length of genome: %lu\n", total_len );
  printf( "# Number of sequences: %d\n", total_records );
  printf( "# Base composition: A C G T N other\n" );
  printf( "%d %d %d %d %d %d\n\n",
	  base_comp[0], base_comp[1],
	  base_comp[2], base_comp[3],
	  base_comp[4], base_comp[5] );

  printf( "# Dinucleotide composition: AA AC AG AT CA CC...\n" );
  for( i = 0; i < 16; i++ ) {
    printf( "%d ", dinuc_comp[i] );
  }
  printf( "\n\n" );
  
  printf( "# Distribution of gap lengths, i.e., N-runs (log10)\n" );
  for( i = 0; i <10; i++ ) {
    printf( "%d ", N_runs[i] );
  }
  printf( "\n\n" );

  qsort( &lengths[0], total_records, sizeof( size_t ), cmpintp );
  i = 0;
  while( n50_tmp < (total_len)/2 ) {
    n50_tmp += lengths[i++];
  }
  printf( "# N50 sequence is number %d and length: %d\n\n",
	  (int)(i-1), (int)lengths[i-1] );
  
  printf( "# 20 longest sequence lengths:\n" );
  for( i = 0; i < 20; i++ ) {
    printf( "%d\n", (int)lengths[i] );
  }

    
}

static int cmpintp( const void *p1, const void *p2 ) {
  int len1, len2;
  len1 = *(int*)p1;
  len2 = *(int*)p2;
  if ( len1 < len2 ) {
    return 1;
  }
  if ( len1 > len2 ) {
    return -1;
  }
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

