#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys

import cyvcf2
import numpy

PERCENTILES=(90, 99, 99.9)
GL_HET=(0.35, 0.65)

def main(tumour, filter_germline_het, pass_only, just_best, info_af, af_name, min_af, min_dp, vcf_filename, out_fh):
  logging.info('reading from %s...', vcf_filename)

  # we just want to get all the AFs
  afs = []
  vcf = cyvcf2.VCF(vcf_filename)
  sample_id = vcf.samples.index(tumour)
  germline_id = 1 if sample_id == 0 else 0
  logging.debug('sample_id %i germline_id %i', sample_id, germline_id)
  skipped = 0
  gaf_range = (1.0, 0.0)

  for v in cyvcf2.VCF('-'):
    if pass_only and v.FILTER is not None:
      logging.debug('skipping non-pass at %s:%s', v.CHROM, v.POS)
      skipped += 1
      continue
    # gl af
    if not info_af:
      gaf = v.format(af_name)[germline_id][0] # mutect2 af
      gaf_range = (min([gaf_range[0], gaf]), max([gaf_range[1], gaf]))
      logging.debug('gaf %s range %s', gaf, gaf_range)
      if filter_germline_het and GL_HET[0] < gaf < GL_HET[1]:
        logging.debug('skipping germline het at %s:%s', v.CHROM, v.POS)
        skipped += 1
        continue

    # tumour af
    if info_af:
      af = v.INFO[af_name] # calculated af
    else:
      af = v.format(af_name)[sample_id] # mutect2 af

    if af < min_af:
      logging.debug('skipping low AF %f at %s:%s', af, v.CHROM, v.POS)
      skipped += 1
      continue

    if min_dp > 0:
      dp = v.format("DP")[sample_id]
      if dp < min_dp:
        logging.debug('skipping low DP %f at %s:%s', dp, v.CHROM, v.POS)
        skipped += 1
        continue

    logging.debug('appending %s to afs', af)
    afs.append(af)

  if len(afs) == 0:
    logging.warn('No afs')
    answer = [0.0, 0.5, 1.0]
  else:
    logging.debug('%i afs: %s', len(afs), afs)
    answer = numpy.percentile(afs, PERCENTILES)

  if out_fh is not None:
    if just_best:
      out_fh.write('{:.2f}'.format(answer[1]))
    else:
      out_fh.write('Lower\tBest\tUpper\n')
      out_fh.write('{:.2f}\t{:.2f}\t{:.2f}\n'.format(answer[0], answer[1], answer[2]))

  logging.info('done. skipped %i included %i. gaf range %s', skipped, len(afs), gaf_range)

  return {'lower': answer[0], 'best': answer[1], 'upper': answer[2]}

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Read a VCF, estimate purity')
  parser.add_argument('--tumour', help='name of tumour')
  parser.add_argument('--info_af', action='store_true', help='af is in info field')
  parser.add_argument('--af_name', default='AF', help='af is in info field')
  parser.add_argument('--pass_only', action='store_true', help='just pass variants')
  parser.add_argument('--just_best', action='store_true', help='only write best estimate')
  parser.add_argument('--min_af', type=float, default=0, required=False, help='only consider variants with af of at least this')
  parser.add_argument('--min_dp', type=float, default=0, required=False, help='only consider variants with dp of at least this')
  parser.add_argument('--filter_germline_het', action='store_true', help='exclude hets in germline (possible LOH)')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.tumour, args.filter_germline_het, args.pass_only, args.just_best, args.info_af, args.af_name, args.min_af, args.min_dp, '-', sys.stdout)
