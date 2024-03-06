
######## Set title ########
fgRed=$(tput setaf 1)
txBold=$(tput bold)
txReset=$(tput sgr0)

termwidth="$(tput cols)"
mssg="Welcome to the MEMDS analysis pipeline"
let padding=($termwidth - ${#mssg} - 3)/2

printf "\n%${padding}s ${fgRed}${txBold} $mssg %${padding}s ${txReset}\n"

######## Gather path infromation ########
echo "Please provide path(s) to the MEMDS pipeline script folder. E.g: My_computer/My_experiment/scripts"
printf "Multiple paths can be provided, separated by semi-colon (;)\n\n"

read path_str

paths=$(echo "$path_str" | tr ";" "\n")

printf "\n${fgRed}${txBold}These paths were provided:${txReset}\n"

for path in $paths; do
    echo "$path"
done

printf "\nAre these paths correct? Y/N: \n"
while read input; do
    case $input in
        [yY])
            echo 'Continuing'
            break
            ;;
        [nN])
            echo "Incorrect paths encountered. Exiting"
            exit
            ;;
        *)
            printf "Please input Y or N" >&2
    esac
done

######## List possible steps ########
explanation="
${fgRed}${txBold}Choose one of the following steps:${txReset}\n
1. Merge partial fastq files (optional)\n
2. Prepare pipeline parameter table\n
3. Quality control of the experimental data\n
4. Trim and validate read barcodes\n
5. Sort reads by origin gene\n
6. Map reads to the reference sequence\n
7. Create files for IGV inspection\n
8. Identify potential mutations\n
9. Identify sequencing quality at mutated positions\n
10. Create consensus mutation tables\n"

printf "$explanation"

######## Execute the steps ########
for path in $paths; do

    ########
    printf "\n^v^v^v^v^v^v^v^v^v^\n"
    printf "${txBold}Working on $path${txReset}\n"
    echo "Please choose step to run [1-10]: "

    read step

    while true; do
        if [[ ! $step =~ ^[0-9]+$ ]] || [ $step -gt 10 ]; then
            echo "Please choose ${txBold}valid${txReset} step to run [1-10]: "
            read step
            continue
        else
            break
        fi
    done

    while true; do
        printf "\nIs the step correct? Y/N: \n"

        read input
        case $input in
            [yY])
                echo "Executing $step"
                break
                ;;
            [nN])
                echo "Please choose step to run [1-10]: "
                read step

                while true; do
                    if [[ ! $step =~ ^[0-9]+$ ]] || [ $step -gt 10 ]; then
                        echo "Please choose ${txBold}valid${txReset} step to run [1-10]: "
                        read step
                        continue
                    else
                        break
                    fi
                done

                continue
                ;;
            *)
                printf "Please input Y or N" >&2
        esac
    done
    ########

    case $step in
        1)
            printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. Create file merging script\n
2. Merge partial fastq files\n\n"
            printf "Please choose sub-step to run [1-2]: "

            read substep

            while true; do
                if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 2 ]; then
                    echo "Please choose ${txBold}valid${txReset} sub-step to run [1-2]: "
                    read substep
                    continue
                else
                    break
                fi
            done

            case $substep in
                1)
                    cd $path
                    bash concatenate_partfiles.sh
                    ;;
                2)
                    cd $path
                    srun bash fastq_merging/samples_table.concat.sh
            esac
            ;;
        2)
            printf "\n${fgRed}${txBold}Choose mode: PE (Paired-end) | SE (single-end)${txReset}\n";

            read mode
            until [ $mode == "PE" ] || [ $mode == "SE" ]; do
                printf "Please choose ${txBold}PE${txReset} or ${txBold}SE${txReset}: \n"
                read mode
            done

            cd $path
            bash setting_1-${mode}.sh
            ;;
        3)
            printf "\n${fgRed}${txBold}Choose mode: PE (Paired-end) | SE (single-end)${txReset}\n";

            read mode
            until [ $mode == "PE" ] || [ $mode == "SE" ]; do
                printf "Please choose ${txBold}PE${txReset} or ${txBold}SE${txReset}: \n"
                read mode
            done

            printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. PE/SE: Assess raw data quality with FastQC\n
2. PE: Merge paired-end reads with PEAR
   SE: Quality-filter raw data with Cutadapt+Trimmomatic\n
3. PE: Quality-filter raw data with Cutadapt+Trimmomatic
   SE: Assess quality-filtered data with FastQC\n
4. PE: Assess quality-filtered data with FastQC\n\n"
            printf "Please choose sub-step to run [1-4]: "

            read substep

            while true; do
                if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 4 ]; then
                    echo "Please choose ${txBold}valid${txReset} sub-step to run [1-4]: "
                    read substep
                    continue
                else
                    break
                fi
            done

            cd $path
            bash filter-${mode}4.sh $substep
            ;;
        4)
            cd $path
            bash trim7.sh 1
            ;;
        5)
            cd $path
            bash sort2.sh 1
            ;;
        6)
            printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. Prepare reference file for alignment\n
2. Align experimental data to the reference\n\n"
            printf "Please choose sub-step to run [1-2]: "

            read substep

            while true; do
                if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 2 ]; then
                    echo "Please choose ${txBold}valid${txReset} sub-step to run [1-2]: "
                    read substep
                    continue
                else
                    break
                fi
            done

            cd $path
            bash bwa9.sh $substep
            ;;
        7)
            cd $path
            bash create_dummy_genome5.sh 1
            ;;
        8)
            cd $path
            bash sam_to_mutation-table_5.2.sh 1
            ;;
        9)
            cd $path
            bash sam_to_mutation-list-3.sh 1
            ;;
        10)
            printf "\n${fgRed}${txBold}Choose one of the following sub-steps:${txReset}\n
1. Summarise mutations per read\n
2. Create consensus mutation tables per defined cut-offs\n\n"
            printf "Please choose sub-step to run [1-2]: "

            read substep

            while true; do
                if [[ ! $substep =~ ^[0-9]+$ ]] || [ $substep -gt 2 ]; then
                    echo "Please choose ${txBold}valid${txReset} sub-step to run [1-2]: "
                    read substep
                    continue
                else
                    break
                fi
            done

            cd $path
            bash consensus_15.1.sh $substep
    esac
done
