#!/usr/bin/env nextflow
// NXF_OFFLINE=true nextflow run main.nf -resume  -with-trace -with-report
nextflow.enable.dsl=2

// === PROCESSUS 1 : Conversion pVCF → PLINK ===
process pVCF_to_plink {
  tag "$vcf_file"

  input:
  path vcf_file

  output:
  tuple path("*.bed"), path("*.bim"), path("*.fam"), emit: bfiles

  script:
  """
  bash ${params.scriptDir}/pVCF_to_plinkNF.sh ${vcf_file} .
  """
}

// === PROCESSUS 2 : Filtres + merge PLINK ===
process plink_filter_merge {
  tag "plink_filter_merge"

  // on utilise le staging natif Nextflow avec path("bfiles/*").
  // Tous les .bed/.bim/.fam sont stagés dans un sous-dossier "bfiles/" du workdir.
  input:
  path("bfiles/*")

  output:
  path "merged_dataset.bed",                   emit: merged_bed
  path "merged_dataset.bim",                   emit: merged_bim
  path "merged_dataset.fam",                   emit: merged_fam
  path "stats/merged_dataset.missing.smiss",   emit: merged_smiss
  path "stats/merged_dataset.missing.vmiss",   emit: merged_vmiss
  path "stats/merged_dataset_unrelated.kin0",  emit: merged_unrelated

  publishDir "${params.outputDir}/results", mode: 'link'

  script:
  """
  bash "${params.scriptDir}/plink_filters_and_mergeNF.sh" "bfiles" . ${params.PLINK_filter_options}
  """
}

process select_unrel_topN {
  tag "select_unrel_topN (topN=${topN})"

  input:
  path kin0_file
  path smiss_file
  val  topN

  output:
  path "keep_unrel_topN.iid.txt", emit: keep_unrel_topN

  publishDir "${params.outputDir}/results", mode: 'link'

  script:
  """
  bash ${params.scriptDir}/select_unrel_topN.sh ${kin0_file} ${smiss_file} ${topN} keep_unrel_topN.iid.txt
  """
}

// === PROCESSUS 4 : gVCF → LRR/BAF  ===
process gvcf_to_signalintensity {
  tag { gvcf_list.simpleName }

  input:
  path gvcf_list
  path merged_bim

  output:
  path "LRR_BAF/*.baf_lrr.tsv", emit: baflrr_files

  // Chaque job publie ses fichiers dans le même répertoire LRR_BAF
  publishDir "${params.outputDir}/results/LRR_BAF", mode: 'link', saveAs: { it.tokenize('/').last() }

  script:
  """
  mkdir -p LRR_BAF
  bash ${params.scriptDir}/multithread_baf_lrr.sh --batch_list ${gvcf_list} --output_dir LRR_BAF --bim ${merged_bim}
  """
}

process split_iids {
  tag "split_iids (n_batches=${n_batches})"

  input:
  path iid_list
  path("lrr_baf_all/*")   // tous les fichiers BAF des batches stagés dans lrr_baf_all/
  val  n_batches

  output:
  path "batches/iid_batch_*", emit: iid_batches

  publishDir "${params.outputDir}/results", mode: 'link'

  script:
  """
  set -euo pipefail
  bash ${params.scriptDir}/make_iid_batches.sh \
       "${iid_list}" \
       "lrr_baf_all" \
       "${n_batches}" \
       .
  """
}

// pfb_partial sur chaque batch (Nextflow parallélise nativement)
process pfb_partial {
  tag { "pfb_partial(${batch.simpleName})" }

  input:
  tuple path(batch), val(baf_dir)

  output:
  path "partial/*.parquet", emit: partials

  publishDir "${params.outputDir}/results", mode: 'link'

  script:
  """
  set -euo pipefail
  mkdir -p partial
  export POLARS_MAX_THREADS=${params.polars_threads}

  python3 ${params.scriptDir}/pfb_partial.py \
    "${batch}" \
    "${baf_dir}" \
    "partial/\$(basename "${batch}").parquet"
  """
}

// Réduction finale : agrège tous les partiels → PFB.tsv
process pfb_reduce {
  tag "pfb_reduce"

  input:
  path partials
  val  out_name

  output:
  path out_name, emit: pfb_tsv

  publishDir "${params.outputDir}/results", mode: 'link'

  script:
  """
  set -euo pipefail
  python3 ${params.scriptDir}/pfb_reduce.py "${out_name}" .
  """
}

process filter_baf_lrr_files {
  // filter
  tag "filter_baf_lrr_batch"

  input:
  path baf_files   // batch de fichiers BAF stagés dans le workdir
  path pfb_file

  output:
  path "*_filtered.baf_lrr.tsv", emit: filtered_baf

  publishDir "${params.outputDir}/results/LRR_BAF_filtered", mode: 'link'

  script:
  """
  set -euo pipefail
  for baf_file in *.baf_lrr.tsv; do
    sample_id="\${baf_file%.baf_lrr.tsv}"
    awk -F'\t' -f ${params.scriptDir}/filter_BAFLRR_by_pfbFinal.awk \
      "${pfb_file}" "\${baf_file}" > "\${sample_id}_filtered.baf_lrr.tsv"
  done
  """
}


// === WORKFLOW ===
workflow {

  main:
    if (!params.inputDirgVCF) {
        error "You must provide a directory with gVCF files"
    }
    if (!params.inputDirpVCF) {
        error "You must provide a directory with pVCF files"
    }

    // 1) pVCF → PLINK (un job par VCF, parallélisé nativement)
    Channel
      .fromPath("${params.inputDirpVCF}/*.vcf.gz", checkIfExists: true)
      .filter { vcf -> !(vcf.name =~ /chr[XY]\.vcf\.gz$/) }
      .set { vcf_files }

    converted = pVCF_to_plink(vcf_files)

    // 2) Collecter tous les .bed/.bim/.fam et les passer à plink_filter_merge.
    //    Le staging natif path("bfiles/*") gère tout dans le workdir du processus.
    converted
      .flatten()
      .collect()
      .set { all_bfiles }

    plink_filter_merge(all_bfiles)

    // 3) Sélection des individus non-apparentés
    select_unrel_topN(
      plink_filter_merge.out.merged_unrelated,
      plink_filter_merge.out.merged_smiss,
      params.topN ?: 700
    )

    // 4) Découper les gVCF en batches → 1 job SLURM par batch de ~gvcf_batch_size individus
    Channel
      .fromPath("${params.inputDirgVCF}/*.gvcf.gz", checkIfExists: true)
      .map { vcf ->
        def sid = vcf.getBaseName().replaceAll(/\.(gvcf|vcf)(\.gz)?$/, "")
        return "${sid}\t${vcf.toAbsolutePath()}"
      }
      .buffer(size: params.gvcf_batch_size ?: 125, remainder: true)
      .map { lines ->
        def uuid = UUID.randomUUID().toString().take(8)
        def f = file("gvcf_batch_${uuid}.tsv")
        f.text = lines.join("\n")
        return f
      }
      .set { gvcf_batch_files }

    // 5)  jobs SLURM en parallèle 
    gvcf_to_signalintensity(gvcf_batch_files, plink_filter_merge.out.merged_bim)

    iid_ch = select_unrel_topN.out.keep_unrel_topN

    // Collecter tous les fichiers BAF de tous les batches
    all_baf_ch = gvcf_to_signalintensity.out.baflrr_files.flatten().collect()

    // 6) split_iids : tous les BAF stagés ensemble dans lrr_baf_all/
    split_iids(
      iid_ch,
      all_baf_ch,
      params.n_batches
    )

    // 7) pfb_partial en parallèle sur chaque batch
    //    On utilise le chemin publié (stable) pour baf_dir afin d'éviter
    //    toute dépendance sur le workdir temporaire de gvcf_to_signalintensity.
    def baf_dir_path = params.baf_dir_path ?: "${params.outputDir}/results/LRR_BAF"

    split_iids.out.iid_batches
      .flatten()
      .map { batch -> tuple(batch, baf_dir_path) }
      | pfb_partial

    // 8) pfb_reduce : agrège tous les partiels → PFB.tsv
    pfb_partial.out.partials
      .collect()
      .set { all_partials }

    pfb_reduce(
      all_partials,
      "PFB.tsv"
    )

    // 9) Filtrer par batches 
    gvcf_to_signalintensity.out.baflrr_files
      .flatten()
      .collate(params.filter_batch_size ?: 250)
      .set { baf_batches_ch }

    filter_baf_lrr_files(baf_batches_ch, pfb_reduce.out.pfb_tsv)
}
