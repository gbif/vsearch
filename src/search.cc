/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.


  GNU General Public License version 3

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.


  The BSD 2-Clause License

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

*/

#include "vsearch.h"
#include "align_simd.h"
#include "dbindex.h"
#include "mask.h"
#include "minheap.h"
#include "otutable.h"
#include "udb.h"
#include "unique.h"
#include "utils/fatal.hpp"
#include "utils/maps.hpp"
#include "utils/xpthread.hpp"
#include <algorithm>  // std::min
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint> // uint64_t, int64_t
#include <cstdio>  // std::FILE, std::fprintf, std::fclose, std::size_t
#include <cstring>  // std::strlen, std::memset, std::strcpy
#include <pthread.h>
#include <vector>


/* db_open tracks whether the database has been loaded; transitions
   false→true exactly once (before the server accepts requests) and is
   never written again, so it is safe as a static. */
static bool db_open = false;

struct SearchContext {
  // Thread management
  struct searchinfo_s *si_plus         = nullptr;
  struct searchinfo_s *si_minus        = nullptr;
  pthread_t           *pthread_arr     = nullptr;
  pthread_attr_t       attr            {};

  // Query input
  fastx_handle         query_fastx_h   = nullptr;

  // Search constants (set in search_prep)
  int tophits  = 0;
  int seqcount = 0;

  // Synchronization
  pthread_mutex_t mutex_input  {};
  pthread_mutex_t mutex_output {};

  // Statistics
  int      qmatches           = 0;
  uint64_t qmatches_abundance = 0;
  int      queries            = 0;
  uint64_t queries_abundance  = 0;

  // DB hit tracking
  uint64_t *dbmatched = nullptr;

  // Output file pointers (16)
  FILE *fp_samout            = nullptr;
  FILE *fp_alnout            = nullptr;
  FILE *fp_userout           = nullptr;
  FILE *fp_blast6out         = nullptr;
  FILE *fp_uc                = nullptr;
  FILE *fp_fastapairs        = nullptr;
  FILE *fp_matched           = nullptr;
  FILE *fp_notmatched        = nullptr;
  FILE *fp_dbmatched         = nullptr;
  FILE *fp_dbnotmatched      = nullptr;
  FILE *fp_otutabout         = nullptr;
  FILE *fp_mothur_shared_out = nullptr;
  FILE *fp_biomout           = nullptr;
  FILE *fp_lcaout            = nullptr;
  FILE *fp_qsegout           = nullptr;
  FILE *fp_tsegout           = nullptr;

  // Output line counters
  int count_matched    = 0;
  int count_notmatched = 0;
};

struct ThreadArgs {
  int64_t       t;
  SearchContext *ctx;
};


auto search_output_results(std::vector<struct hit> const & hits,
                           char const * query_head,
                           int qseqlen,
                           char const * qsequence,
                           char const * qsequence_rc,
                           int qsize,
                           SearchContext *ctx) -> void
{
  xpthread_mutex_lock(&ctx->mutex_output);

  /* show results */
  auto const toreport = std::min(opt_maxhits, static_cast<int64_t>(hits.size()));

  if (ctx->fp_alnout != nullptr)
    {
      results_show_alnout(ctx->fp_alnout,
                          hits.data(),
                          toreport,
                          query_head,
                          qsequence,
                          qseqlen);
    }

  if (ctx->fp_lcaout != nullptr)
    {
      results_show_lcaout(ctx->fp_lcaout,
                          hits.data(),
                          toreport,
                          query_head);
    }

  if (ctx->fp_samout != nullptr)
    {
      results_show_samout(ctx->fp_samout,
                          hits.data(),
                          toreport,
                          query_head,
                          qsequence,
                          qsequence_rc);
    }

  if (toreport != 0)
    {
      double const top_hit_id = hits[0].id;

      if ((opt_otutabout != nullptr) || (opt_mothur_shared_out != nullptr) || (opt_biomout != nullptr))
        {
          otutable_add(query_head,
                       db_getheader(hits[0].target),
                       qsize);
        }

      for (auto t = 0; t < toreport; t++)
        {
          auto const * hp = &hits[t];

          if ((opt_top_hits_only != 0) && (hp->id < top_hit_id))
            {
              break;
            }

          if (ctx->fp_fastapairs != nullptr)
            {
              results_show_fastapairs_one(ctx->fp_fastapairs,
                                          hp,
                                          query_head,
                                          qsequence,
                                          qsequence_rc);
            }

          if (ctx->fp_qsegout != nullptr)
            {
              results_show_qsegout_one(ctx->fp_qsegout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (ctx->fp_tsegout != nullptr)
            {
              results_show_tsegout_one(ctx->fp_tsegout,
                                       hp);
            }

          if (ctx->fp_uc != nullptr)
            {
              if ((t==0) || (opt_uc_allhits != 0))
                {
                  results_show_uc_one(ctx->fp_uc,
                                      hp,
                                      query_head,
                                      qseqlen,
                                      hp->target);
                }
            }

          if (ctx->fp_userout != nullptr)
            {
              results_show_userout_one(ctx->fp_userout,
                                       hp,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (ctx->fp_blast6out != nullptr)
            {
              results_show_blast6out_one(ctx->fp_blast6out,
                                         hp,
                                         query_head,
                                         qseqlen);
            }
        }
    }
  else
    {
      if ((opt_otutabout != nullptr) || (opt_mothur_shared_out != nullptr) || (opt_biomout != nullptr))
        {
          otutable_add(query_head,
                       nullptr,
                       qsize);
        }

      if (ctx->fp_uc != nullptr)
        {
          results_show_uc_one(ctx->fp_uc,
                              nullptr,
                              query_head,
                              qseqlen,
                              0);
        }

      if (opt_output_no_hits != 0)
        {
          if (ctx->fp_userout != nullptr)
            {
              results_show_userout_one(ctx->fp_userout,
                                       nullptr,
                                       query_head,
                                       qsequence,
                                       qseqlen,
                                       qsequence_rc);
            }

          if (ctx->fp_blast6out != nullptr)
            {
              results_show_blast6out_one(ctx->fp_blast6out,
                                         nullptr,
                                         query_head,
                                         qseqlen);
            }
        }
    }

  if (not hits.empty())
    {
      ctx->count_matched++;
      if (opt_matched != nullptr)
        {
          fasta_print_general(ctx->fp_matched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              strlen(query_head),
                              qsize,
                              ctx->count_matched,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }
    }
  else
    {
      ctx->count_notmatched++;
      if (opt_notmatched != nullptr)
        {
          fasta_print_general(ctx->fp_notmatched,
                              nullptr,
                              qsequence,
                              qseqlen,
                              query_head,
                              strlen(query_head),
                              qsize,
                              ctx->count_notmatched,
                              -1.0,
                              -1, -1, nullptr, 0.0);
        }
    }

  /* update matching db sequences */
  for (auto const & hit : hits) {
    if (hit.accepted or hit.weak) {
      ctx->dbmatched[hit.target] += opt_sizein ? qsize : 1;
    }
  }

  xpthread_mutex_unlock(&ctx->mutex_output);
}


auto search_query(int64_t t, SearchContext *ctx) -> int
{
  for (int s = 0; s < opt_strand; s++)
    {
      struct searchinfo_s * si = (s != 0) ? ctx->si_minus + t : ctx->si_plus + t;

      /* mask query */
      if (opt_qmask == MASK_DUST)
        {
          dust(si->qsequence, si->qseqlen);
        }
      else if ((opt_qmask == MASK_SOFT) && (opt_hardmask != 0))
        {
          hardmask(si->qsequence, si->qseqlen);
        }

      /* perform search */
      search_onequery(si, opt_qmask);
    }

  std::vector<struct hit> hits;

  search_joinhits(ctx->si_plus + t,
                  opt_strand > 1 ? ctx->si_minus + t : nullptr,
                  hits);

  search_output_results(hits,
                        ctx->si_plus[t].query_head,
                        ctx->si_plus[t].qseqlen,
                        ctx->si_plus[t].qsequence,
                        opt_strand > 1 ? ctx->si_minus[t].qsequence : nullptr,
                        ctx->si_plus[t].qsize,
                        ctx);

  /* free memory for alignment strings */
  for (auto const & hit : hits) {
    if (hit.aligned) {
      xfree(hit.nwalignment);
    }
  }

  return static_cast<int>(hits.size());
}


auto search_thread_run(int64_t t, SearchContext *ctx) -> void
{
  while (true)
    {
      xpthread_mutex_lock(&ctx->mutex_input);

      if (fastx_next(ctx->query_fastx_h,
                     (opt_notrunclabels == 0),
                     chrmap_no_change_vector.data()))
        {
          char const * qhead = fastx_get_header(ctx->query_fastx_h);
          int const query_head_len = fastx_get_header_length(ctx->query_fastx_h);
          char const * qseq = fastx_get_sequence(ctx->query_fastx_h);
          int const qseqlen = fastx_get_sequence_length(ctx->query_fastx_h);
          int const query_no = fastx_get_seqno(ctx->query_fastx_h);
          int const qsize = fastx_get_abundance(ctx->query_fastx_h);

          for (int s = 0; s < opt_strand; s++)
            {
              struct searchinfo_s * si = (s != 0) ? ctx->si_minus + t : ctx->si_plus + t;

              si->query_head_len = query_head_len;
              si->qseqlen = qseqlen;
              si->query_no = query_no;
              si->qsize = qsize;
              si->strand = s;

              /* allocate more memory for header and sequence, if necessary */

              if (si->query_head_len + 1 > si->query_head_alloc)
                {
                  si->query_head_alloc = si->query_head_len + 2001;
                  si->query_head = (char *)
                    xrealloc(si->query_head, (size_t) (si->query_head_alloc));
                }

              if (si->qseqlen + 1 > si->seq_alloc)
                {
                  si->seq_alloc = si->qseqlen + 2001;
                  si->qsequence = (char *)
                    xrealloc(si->qsequence, (size_t) (si->seq_alloc));
                }
            }

          /* plus strand: copy header and sequence */
          strcpy(ctx->si_plus[t].query_head, qhead);
          strcpy(ctx->si_plus[t].qsequence, qseq);

          /* get progress as amount of input file read */
          uint64_t const progress = fastx_get_position(ctx->query_fastx_h);

          /* let other threads read input */
          xpthread_mutex_unlock(&ctx->mutex_input);

          /* minus strand: copy header and reverse complementary sequence */
          if (opt_strand > 1)
            {
              strcpy(ctx->si_minus[t].query_head, ctx->si_plus[t].query_head);
              reverse_complement(ctx->si_minus[t].qsequence,
                                 ctx->si_plus[t].qsequence,
                                 ctx->si_plus[t].qseqlen);
            }

          int const match = search_query(t, ctx);

          /* lock mutex for update of global data and output */
          xpthread_mutex_lock(&ctx->mutex_output);

          /* update stats */
          ++ctx->queries;
          ctx->queries_abundance += qsize;

          if (match != 0)
            {
              ++ctx->qmatches;
              ctx->qmatches_abundance += qsize;
            }

          /* show progress */
          progress_update(progress);

          xpthread_mutex_unlock(&ctx->mutex_output);
        }
      else
        {
          xpthread_mutex_unlock(&ctx->mutex_input);
          break;
        }
    }
}


auto search_thread_init(struct searchinfo_s * si, SearchContext *ctx) -> void
{
  /* thread specific initialiation */
  si->uh = unique_init();
  si->kmers = (count_t *) xmalloc((ctx->seqcount * sizeof(count_t)) + 32);
  si->m = minheap_init(ctx->tophits);
  si->hits = (struct hit *) xmalloc
    (sizeof(struct hit) * (ctx->tophits) * opt_strand);
  si->qsize = 1;
  si->query_head_alloc = 0;
  si->query_head = nullptr;
  si->seq_alloc = 0;
  si->qsequence = nullptr;
  si->s = search16_init(opt_match,
                        opt_mismatch,
                        opt_gap_open_query_left,
                        opt_gap_open_target_left,
                        opt_gap_open_query_interior,
                        opt_gap_open_target_interior,
                        opt_gap_open_query_right,
                        opt_gap_open_target_right,
                        opt_gap_extension_query_left,
                        opt_gap_extension_target_left,
                        opt_gap_extension_query_interior,
                        opt_gap_extension_target_interior,
                        opt_gap_extension_query_right,
                        opt_gap_extension_target_right);
}


auto search_thread_exit(struct searchinfo_s * si) -> void
{
  /* thread specific clean up */
  search16_exit(si->s);
  unique_exit(si->uh);
  xfree(si->hits);
  minheap_exit(si->m);
  xfree(si->kmers);
  if (si->query_head != nullptr)
    {
      xfree(si->query_head);
    }
  if (si->qsequence != nullptr)
    {
      xfree(si->qsequence);
    }
}


auto search_thread_worker(void * vp) -> void *
{
  auto *args = static_cast<ThreadArgs *>(vp);
  search_thread_run(args->t, args->ctx);
  return nullptr;
}


auto search_thread_worker_run(SearchContext *ctx) -> void
{
  /* initialize threads, start them, join them and return */

  xpthread_attr_init(&ctx->attr);
  xpthread_attr_setdetachstate(&ctx->attr, PTHREAD_CREATE_JOINABLE);

  auto *thread_args = new ThreadArgs[opt_threads];

  /* init and create worker threads, put them into stand-by mode */
  for (int t = 0; t < opt_threads; t++)
    {
      search_thread_init(ctx->si_plus + t, ctx);
      if (ctx->si_minus != nullptr)
        {
          search_thread_init(ctx->si_minus + t, ctx);
        }
      thread_args[t] = {t, ctx};
      xpthread_create(ctx->pthread_arr + t, &ctx->attr,
                      search_thread_worker, &thread_args[t]);
    }

  /* finish and clean up worker threads */
  for (int t = 0; t < opt_threads; t++)
    {
      xpthread_join(ctx->pthread_arr[t], nullptr);
      search_thread_exit(ctx->si_plus + t);
      if (ctx->si_minus != nullptr)
        {
          search_thread_exit(ctx->si_minus + t);
        }
    }

  delete[] thread_args;
  xpthread_attr_destroy(&ctx->attr);
}


auto search_prep(char * cmdline, char * progheader, SearchContext *ctx) -> void
{
  /* open output files */

  if (opt_alnout != nullptr)
    {
      ctx->fp_alnout = fopen_output(opt_alnout);
      if (ctx->fp_alnout == nullptr)
        {
          fatal("Unable to open alignment output file for writing");
        }

      fprintf(ctx->fp_alnout, "%s\n", cmdline);
      fprintf(ctx->fp_alnout, "%s\n", progheader);
    }

  if (opt_lcaout != nullptr)
    {
      ctx->fp_lcaout = fopen_output(opt_lcaout);
      if (ctx->fp_lcaout == nullptr)
        {
          fatal("Unable to open lca output file for writing");
        }
    }

  if (opt_samout != nullptr)
    {
      ctx->fp_samout = fopen_output(opt_samout);
      if (ctx->fp_samout == nullptr)
        {
          fatal("Unable to open SAM output file for writing");
        }
    }

  if (opt_userout != nullptr)
    {
      ctx->fp_userout = fopen_output(opt_userout);
      if (ctx->fp_userout == nullptr)
        {
          fatal("Unable to open user-defined output file for writing");
        }
    }

  if (opt_blast6out != nullptr)
    {
      ctx->fp_blast6out = fopen_output(opt_blast6out);
      if (ctx->fp_blast6out == nullptr)
        {
          fatal("Unable to open blast6-like output file for writing");
        }
    }

  if (opt_uc != nullptr)
    {
      ctx->fp_uc = fopen_output(opt_uc);
      if (ctx->fp_uc == nullptr)
        {
          fatal("Unable to open uc output file for writing");
        }
    }

  if (opt_fastapairs != nullptr)
    {
      ctx->fp_fastapairs = fopen_output(opt_fastapairs);
      if (ctx->fp_fastapairs == nullptr)
        {
          fatal("Unable to open fastapairs output file for writing");
        }
    }

  if (opt_qsegout != nullptr)
    {
      ctx->fp_qsegout = fopen_output(opt_qsegout);
      if (ctx->fp_qsegout == nullptr)
        {
          fatal("Unable to open qsegout output file for writing");
        }
    }

  if (opt_tsegout != nullptr)
    {
      ctx->fp_tsegout = fopen_output(opt_tsegout);
      if (ctx->fp_tsegout == nullptr)
        {
          fatal("Unable to open tsegout output file for writing");
        }
    }

  if (opt_matched != nullptr)
    {
      ctx->fp_matched = fopen_output(opt_matched);
      if (ctx->fp_matched == nullptr)
        {
          fatal("Unable to open matched output file for writing");
        }
    }

  if (opt_notmatched != nullptr)
    {
      ctx->fp_notmatched = fopen_output(opt_notmatched);
      if (ctx->fp_notmatched == nullptr)
        {
          fatal("Unable to open notmatched output file for writing");
        }
    }

  if (opt_otutabout != nullptr)
    {
      ctx->fp_otutabout = fopen_output(opt_otutabout);
      if (ctx->fp_otutabout == nullptr)
        {
          fatal("Unable to open OTU table (text format) output file for writing");
        }
    }

  if (opt_mothur_shared_out != nullptr)
    {
      ctx->fp_mothur_shared_out = fopen_output(opt_mothur_shared_out);
      if (ctx->fp_mothur_shared_out == nullptr)
        {
          fatal("Unable to open OTU table (mothur format) output file for writing");
        }
    }

  if (opt_biomout != nullptr)
    {
      ctx->fp_biomout = fopen_output(opt_biomout);
      if (ctx->fp_biomout == nullptr)
        {
          fatal("Unable to open OTU table (biom 1.0 format) output file for writing");
        }
    }

  /* check if it may be an UDB file */

  bool const is_udb = udb_detect_isudb(opt_db);

  if (is_udb)
    {
      if (!db_open)  // TODO: this is a hack
        {
          udb_read(opt_db, true, true);
          results_show_samheader(ctx->fp_samout, cmdline, opt_db);
          show_rusage();
          db_open = true;
        }
      ctx->seqcount = db_getsequencecount();
    }
  else
    {
      db_read(opt_db, 0);
      results_show_samheader(ctx->fp_samout, cmdline, opt_db);
      if (opt_dbmask == MASK_DUST)
        {
          dust_all();
        }
      else if ((opt_dbmask == MASK_SOFT) && (opt_hardmask != 0))
        {
          hardmask_all();
        }
      show_rusage();
      ctx->seqcount = db_getsequencecount();
      dbindex_prepare(1, opt_dbmask);
      dbindex_addallsequences(opt_dbmask);
    }

  /* tophits = the maximum number of hits we need to store */

  if ((opt_maxrejects == 0) || (opt_maxrejects > ctx->seqcount))
    {
      opt_maxrejects = ctx->seqcount;
    }

  if ((opt_maxaccepts == 0) || (opt_maxaccepts > ctx->seqcount))
    {
      opt_maxaccepts = ctx->seqcount;
    }

  ctx->tophits = opt_maxrejects + opt_maxaccepts + MAXDELAYED;

  ctx->tophits = std::min(ctx->tophits, ctx->seqcount);
}


auto search_done(bool skipCloseDB, SearchContext *ctx) -> void
{
  /* clean up, global */
  if (!skipCloseDB) {
    dbindex_free();
    db_free();
  }

  if (opt_lcaout != nullptr)
    {
      fclose(ctx->fp_lcaout);
    }
  if (opt_matched != nullptr)
    {
      fclose(ctx->fp_matched);
    }
  if (opt_notmatched != nullptr)
    {
      fclose(ctx->fp_notmatched);
    }
  if (opt_fastapairs != nullptr)
    {
      fclose(ctx->fp_fastapairs);
    }
  if (opt_qsegout != nullptr)
    {
      fclose(ctx->fp_qsegout);
    }
  if (opt_tsegout != nullptr)
    {
      fclose(ctx->fp_tsegout);
    }
  if (ctx->fp_uc != nullptr)
    {
      fclose(ctx->fp_uc);
    }
  if (ctx->fp_blast6out != nullptr)
    {
      fclose(ctx->fp_blast6out);
    }
  if (ctx->fp_userout != nullptr)
    {
      fclose(ctx->fp_userout);
      clean_up(); // free userfields allocation
    }
  if (ctx->fp_alnout != nullptr)
    {
      fclose(ctx->fp_alnout);
    }
  if (ctx->fp_samout != nullptr)
    {
      fclose(ctx->fp_samout);
    }
  show_rusage();
}


auto usearch_global(struct Parameters const & parameters, char * cmdline, char * progheader, char * fastx, bool skipCloseDB) -> void
{
  SearchContext ctx{};

  search_prep(cmdline, progheader, &ctx);

  if (opt_dbmatched != nullptr)
    {
      ctx.fp_dbmatched = fopen_output(opt_dbmatched);
      if (ctx.fp_dbmatched == nullptr)
        {
          fatal("Unable to open dbmatched output file for writing");
        }
    }

  if (opt_dbnotmatched != nullptr)
    {
      ctx.fp_dbnotmatched = fopen_output(opt_dbnotmatched);
      if (ctx.fp_dbnotmatched == nullptr)
        {
          fatal("Unable to open dbnotmatched output file for writing");
        }
    }

  ctx.dbmatched = (uint64_t *) xmalloc(ctx.seqcount * sizeof(uint64_t *));
  std::memset(ctx.dbmatched, 0, ctx.seqcount * sizeof(uint64_t *));

  otutable_init();

  /* prepare reading of queries */
  ctx.qmatches = 0;
  ctx.qmatches_abundance = 0;
  ctx.queries = 0;
  ctx.queries_abundance = 0;
  // Modification for server mode to avoid segment fault
  ctx.query_fastx_h = fastx_open(fastx);

  /* allocate memory for thread info */
  ctx.si_plus = (struct searchinfo_s *) xmalloc(opt_threads *
                                            sizeof(struct searchinfo_s));
  if (opt_strand > 1)
    {
      ctx.si_minus = (struct searchinfo_s *) xmalloc(opt_threads *
                                                 sizeof(struct searchinfo_s));
    }
  else
    {
      ctx.si_minus = nullptr;
    }

  ctx.pthread_arr = (pthread_t *) xmalloc(opt_threads * sizeof(pthread_t));

  /* init mutexes for input and output */
  xpthread_mutex_init(&ctx.mutex_input, nullptr);
  xpthread_mutex_init(&ctx.mutex_output, nullptr);

  progress_init("Searching", fastx_get_size(ctx.query_fastx_h));
  search_thread_worker_run(&ctx);
  progress_done();

  xpthread_mutex_destroy(&ctx.mutex_output);
  xpthread_mutex_destroy(&ctx.mutex_input);

  xfree(ctx.pthread_arr);
  xfree(ctx.si_plus);
  if (ctx.si_minus != nullptr)
    {
      xfree(ctx.si_minus);
    }

  fastx_close(ctx.query_fastx_h);

  if (! opt_quiet)
    {
      fprintf(stderr, "Matching unique query sequences: %d of %d",
              ctx.qmatches, ctx.queries);
      if (ctx.queries > 0)
        {
          fprintf(stderr, " (%.2f%%)", 100.0 * ctx.qmatches / ctx.queries);
        }
      fprintf(stderr, "\n");
      if (opt_sizein)
        {
          fprintf(stderr, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  ctx.qmatches_abundance, ctx.queries_abundance);
          if (ctx.queries_abundance > 0)
            {
              fprintf(stderr, " (%.2f%%)",
                      100.0 * ctx.qmatches_abundance / ctx.queries_abundance);
            }
          fprintf(stderr, "\n");
        }
    }

  if (opt_log != nullptr)
    {
      fprintf(fp_log, "Matching unique query sequences: %d of %d",
              ctx.qmatches, ctx.queries);
      if (ctx.queries > 0)
        {
          fprintf(fp_log, " (%.2f%%)", 100.0 * ctx.qmatches / ctx.queries);
        }
      fprintf(fp_log, "\n");
      if (opt_sizein)
        {
          fprintf(fp_log, "Matching total query sequences: %" PRIu64 " of %"
                  PRIu64,
                  ctx.qmatches_abundance, ctx.queries_abundance);
          if (ctx.queries_abundance > 0)
            {
              fprintf(fp_log, " (%.2f%%)",
                      100.0 * ctx.qmatches_abundance / ctx.queries_abundance);
            }
          fprintf(fp_log, "\n");
        }
    }


  // Add OTUs with no matches to OTU table
  if ((opt_otutabout != nullptr) || (opt_mothur_shared_out != nullptr) || (opt_biomout != nullptr)) {
    for (int64_t i = 0; i < ctx.seqcount; i++) {
      if (ctx.dbmatched[i] == 0U) {
        otutable_add(nullptr, db_getheader(i), 0);
      }
    }
  }

  if (opt_biomout != nullptr)
    {
      otutable_print_biomout(ctx.fp_biomout);
      fclose(ctx.fp_biomout);
    }

  if (opt_otutabout != nullptr)
    {
      otutable_print_otutabout(ctx.fp_otutabout);
      fclose(ctx.fp_otutabout);
    }

  if (opt_mothur_shared_out != nullptr)
    {
      otutable_print_mothur_shared_out(ctx.fp_mothur_shared_out);
      fclose(ctx.fp_mothur_shared_out);
    }

  otutable_done();

  int count_dbmatched = 0;
  int count_dbnotmatched = 0;

  if ((opt_dbmatched != nullptr) || (opt_dbnotmatched != nullptr))
    {
      for (int64_t i = 0; i < ctx.seqcount; i++)
        {
          if (ctx.dbmatched[i] != 0U)
            {
              count_dbmatched++;
              if (opt_dbmatched != nullptr)
                {
                  fasta_print_general(ctx.fp_dbmatched,
                                      nullptr,
                                      db_getsequence(i),
                                      db_getsequencelen(i),
                                      db_getheader(i),
                                      db_getheaderlen(i),
                                      ctx.dbmatched[i],
                                      count_dbmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0);
                }
            }
          else
            {
              count_dbnotmatched++;
              if (opt_dbnotmatched != nullptr)
                {
                  fasta_print_general(ctx.fp_dbnotmatched,
                                      nullptr,
                                      db_getsequence(i),
                                      db_getsequencelen(i),
                                      db_getheader(i),
                                      db_getheaderlen(i),
                                      db_getabundance(i),
                                      count_dbnotmatched,
                                      -1.0,
                                      -1, -1, nullptr, 0.0);
                }
            }
        }
    }

  xfree(ctx.dbmatched);

  if (opt_dbmatched != nullptr)
    {
      fclose(ctx.fp_dbmatched);
    }
  if (opt_dbnotmatched != nullptr)
    {
      fclose(ctx.fp_dbnotmatched);
    }

  search_done(skipCloseDB, &ctx);
}

auto usearch_global(struct Parameters const & parameters, char * cmdline, char * progheader) -> void
{
  usearch_global(parameters, cmdline, progheader, parameters.opt_usearch_global, false); // original behaviour
}

auto usearch_global_server(struct Parameters const & parameters, char * cmdline, char * progheader, char * query_file) -> void
{
  usearch_global(parameters, cmdline, progheader, query_file, true);   // skips the udb opening
}
