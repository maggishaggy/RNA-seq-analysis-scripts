#!/bin/sh

main(){
    ROOT_PATH=/maui/malvika/Projects/2013-03-09-Ana_Eulalio-Human_and_Salmonella_Typhymurium
    CURRENT_PATH=ID-001066-ID-001077/2015-07-29-pipeline_for_miR_deseq_comparison
    create_folders
    create_sub_folders
    get_deseq_file
    select_regulated_deseq_data
    regulated_genes_arm_preference
    create_mir_clusters_by_seed
    seed_based_clustering
    create_excel_files
}

create_folders(){
    for MAIN_PATH in bin input output
    do
        if ! [ -d $MAIN_PATH ]
        then
            mkdir -p $MAIN_PATH
        fi
    done
}

create_sub_folders(){
    for MAIN_PATH in input/deseq_input input/mature_mir_fasta \
        output/deseq_selected_regulated_genes output/seed_based_clustering \
        output/arm_preference
    do
        if ! [ -d $MAIN_PATH ]
        then
            mkdir -p $MAIN_PATH
        fi
    done
}

get_deseq_file(){
    deseq_path=ID-001066-ID-001077/2014-10-02-READemption_RNA_mapping_miRBase_readlength_14-ID-001066-ID-001077/READemption_analysis/output/deseq
    filename1=deseq_Non_infected_vs_GFP_neg
    cp $ROOT_PATH/$deseq_path/$filename1/*_extended.csv input/deseq_input

    filename2=deseq_Non_infected_vs_GFP_pos
    cp $ROOT_PATH/$deseq_path/$filename2/*_extended.csv input/deseq_input
    
    filename3=deseq_Non_infected_vs_WT
    cp $ROOT_PATH/$deseq_path/$filename3/*_extended.csv input/deseq_input
    
    filename4=deseq_GFP_neg_vs_GFP_pos
    cp $ROOT_PATH/$deseq_path/$filename4/*_extended.csv input/deseq_input
}

select_regulated_deseq_data(){
    inpath=input/deseq_input
    for files in $(ls $inpath)
    do
        threshold=0.99
        outpath=output/deseq_selected_regulated_genes
        python bin/select_regulated_deseq_data.py $inpath/$files $threshold $outpath
    done
}

regulated_genes_arm_preference(){
    inpath=output/deseq_selected_regulated_genes
    for files in $(ls $inpath)
    do
        outpath=output/arm_preference
        python bin/arm_preference_of_regulated_genes.py $inpath/$files $outpath
    done
}

create_mir_clusters_by_seed(){
    inpath=input/mature_mir_fasta
    infile=input/mature_mir_fasta/mature.fa
    if ! [ -f $infile ]
    then
        wget -cP $inpath ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
        gunzip $inpath/mature.fa.gz
    fi
    outfile=input/mature_mir_fasta/miR_clusters_by_seed.txt
    python bin/create_mir_cluster_by_seed.py $infile $outfile
}

seed_based_clustering(){
    inpath=output/deseq_selected_regulated_genes
    cluster_ref_file=input/mature_mir_fasta/miR_clusters_by_seed.txt
    outpath=output/seed_based_clustering
    for files in $(ls $inpath)
    do
        python bin/seed_based_clustering_of_regulated_genes.py $cluster_ref_file $inpath/$files $outpath
    done   
}

create_excel_files(){
    SEND_FOLDER=2015-07-28-Charlotte_Clip_vs_transcriptomics
    mkdir $SEND_FOLDER
    cp -r output/* $SEND_FOLDER
    find ${SEND_FOLDER} -name "*csv" -print0 | xargs -0 -n1 -P24 python ~/bin/csv2xlsx.py
    find ${SEND_FOLDER} -name "*csv" -exec rm {} \;
    zip -r ${SEND_FOLDER}.zip ${SEND_FOLDER}
}
main
