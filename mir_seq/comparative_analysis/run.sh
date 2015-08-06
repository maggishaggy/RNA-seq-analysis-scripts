#!/bin/sh

main(){
    seed_threshold=0.99
    logfc_threshold=0.99
    create_folders
    create_sub_folders
    create_species_folder
    #create_empty_files
    #select_deseq_data
    #####compare miRNAs in two analysis if anything gets affected
    common_microRNA_from_two_analysis 
    pairwise_comparison
    three_samples_comparison
    four_samples_comparison
    get_seed_based_analysis_data
    seed_based_cluster_two_samples_comparison
    seed_based_cluster_three_samples_comparison
    seed_based_cluster_four_samples_comparison
    compile_comparison_data
    convert_csv_to_excel
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
    for MAIN_PATH in input/full_transcriptome \
        input/microRNA_only \
        output/full_transcriptome \
        output/microRNA_only \
        output/full_transcriptome_vs_only_micro \
        input/seed_based_clustering_data \
        output/pairwise_comparison \
        output/three_samples_comparison \
        output/four_samples_comparison \
        output/compiled_data
    do
        if ! [ -d $MAIN_PATH ]
        then
            mkdir -p $MAIN_PATH
        fi
    done
}

create_species_folder(){
    for subfolder in full_transcriptome microRNA_only seed_based_clustering_data
    do
        for species in salmonella shigella listeria staphylococcus
        do
            if ! [ -d input/$subfolder/$species ]
            then
                mkdir -p input/$subfolder/$species
            fi
        done
    done
}

create_empty_files(){
    linc_path=/maui/malvika/Projects/2013-03-09-Ana_Eulalio-Human_and_Salmonella_Typhymurium/lincRNA_comparative_analysis/input/all_annotation
    mir_path=/maui/malvika/Projects/2013-03-09-Ana_Eulalio-Human_and_Salmonella_Typhymurium/microRNA_comparative_analysis/input
    for species in salmonella shigella listeria staphylococcus
    do
        for files in $(ls $linc_path/$species)
        do
            for subfolder in full_transcriptome microRNA_only
            do
                if ! [ -f $mir_path/$subfolder/$species/$files ]
                then
                    touch $mir_path/$subfolder/$species/$files
                fi
            done
        done
    done
}

select_deseq_data(){
    for subpath in full_transcriptome microRNA_only
    do
        for folders in $(ls input/$subpath)
        do
            if ! [ -d output/$subpath/$folders ]
            then
                mkdir -p output/$subpath/$folders
            fi
            inpath=input/$subpath/$folders
            outpath=output/$subpath/$folders
            python bin/select_deseq_data.py $inpath $logfc_threshold $outpath
        done
    done
}

common_microRNA_from_two_analysis(){
    for folders in $(ls output/full_transcriptome)
    do
        full_transcriptome_data=output/full_transcriptome/$folders
        mir_data=output/microRNA_only/$folders
        outpath=output/full_transcriptome_vs_only_micro
        if ! [ -d $outpath/$folders ]
        then
            mkdir -p $outpath/$folders
        fi
        python bin/compare_full_transcriptome_vs_only_mirna_annotation.py $full_transcriptome_data $mir_data $outpath/$folders
    done 
}

pairwise_comparison(){
    for subpath in full_transcriptome microRNA_only
    do
        inpath=output/$subpath
        outpath=output/pairwise_comparison/$subpath
        if ! [ -d $outpath ]
        then
            mkdir -p $outpath
        fi
        python bin/microRNA_comparative_analysis.py $inpath $outpath 2
    done
}

three_samples_comparison(){
    for subpath in full_transcriptome microRNA_only
    do
        inpath=output/$subpath
        outpath=output/three_samples_comparison/$subpath
        if ! [ -d $outpath ]
        then
            mkdir -p $outpath
        fi
        python bin/microRNA_comparative_analysis.py $inpath $outpath 3
    done
}
    
four_samples_comparison(){
    for subpath in full_transcriptome microRNA_only
    do
        inpath=output/$subpath
        outpath=output/four_samples_comparison/$subpath
        if ! [ -d $outpath ]
        then
            mkdir -p $outpath
        fi
        python bin/microRNA_comparative_analysis.py $inpath $outpath 4
    done
}

get_seed_based_analysis_data(){
    main_path=/maui/malvika/Projects/2013-03-09-Ana_Eulalio-Human_and_Salmonella_Typhymurium
    subpath=2015-07-29-pipeline_for_miR_deseq_comparison/output/seed_based_clustering
    data_path=input/seed_based_clustering_data
    sal_path=$main_path/ID-001066-ID-001077/$subpath
    shig_path=$main_path/ID-001752-ID-001759/$subpath
    list_path=$main_path/ID-002411-ID-002418/$subpath
    staph_path=$main_path/ID-003473-ID-003478/$subpath
    cp $sal_path/* $data_path/salmonella
    cp $shig_path/* $data_path/shigella
    cp $list_path/* $data_path/listeria
    cp $staph_path/* $data_path/staphylococcus
}

seed_based_cluster_two_samples_comparison(){
    inpath=input/seed_based_clustering_data
    outpath=output/pairwise_comparison/seed_based_clustering_comparison
    if ! [ -d $outpath ]
    then
        mkdir -p $outpath
    fi
    python bin/seed_based_cluster_comparison.py $inpath $outpath 2 $seed_threshold
}

seed_based_cluster_three_samples_comparison(){
    inpath=input/seed_based_clustering_data
    outpath=output/three_samples_comparison/seed_based_clustering_comparison
    if ! [ -d $outpath ]
    then
        mkdir -p $outpath
    fi
    python bin/seed_based_cluster_comparison.py $inpath $outpath 3 $seed_threshold
}

seed_based_cluster_four_samples_comparison(){
    inpath=input/seed_based_clustering_data
    outpath=output/four_samples_comparison/seed_based_clustering_comparison
    if ! [ -d $outpath ]
    then
        mkdir -p $outpath
    fi
    python bin/seed_based_cluster_comparison.py $inpath $outpath 4 $seed_threshold
}

compile_comparison_data(){
    for subpath in seed_based_clustering_comparison #full_transcriptome microRNA_only 
    do
        outpath=output/compiled_data/$subpath
        if ! [ -d $outpath ]
        then
            mkdir -p $outpath
        fi
        for comp in pairwise_comparison three_samples_comparison four_samples_comparison
        do
            inpath=output/$comp/$subpath
            outfile=$outpath'/compiled_'$comp'.csv'
            rm $outfile
            touch $outfile
            for files in $(ls $inpath)
            do
                echo $files >> $outfile
                cat $inpath/$files >> $outfile
                echo '\n\n' >> $outfile
            done
        done
    done
}

convert_csv_to_excel(){
    SEND_FOLDER=2015-07-31-Ana_Eulalio_microRNA_comparative_study
    mkdir $SEND_FOLDER
    cp -r output/compiled_data/* $SEND_FOLDER
    find ${SEND_FOLDER} -name "*csv" -print0 | xargs -0 -n1 -P24 python ~/bin/csv2xlsx.py
    find ${SEND_FOLDER} -name "*csv" -exec rm {} \;
    zip -r ${SEND_FOLDER}.zip ${SEND_FOLDER}
}

main
