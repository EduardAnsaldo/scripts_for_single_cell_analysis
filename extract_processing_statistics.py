#!/usr/bin/env python
import sys

file_dict = {}
output_file = open('post_aggregation_metric_summary.txt', 'w')

number_of_lanes = sys.argv[1]

for n in range(1,int(number_of_lanes)):
    lane = 'GEM' + str(n)
    file = "./" + lane + "/outs/per_sample_outs/" + lane + "/metrics_summary.csv"
    file_dict[lane] = file

RNA_Q30_percentages = []
barcode_UMI_Q30_percentages = []
median_genes_per_cell = []
confidently_mapped_to_transcriptome = []
mean_antibody_reads_usable_per_cell = []
gene_expression_sequencing_saturation = []
VDJ_mean_used_reads_per_cell = []
gene_expression_number_of_reads = 0
VDJ_number_of_reads = 0
antibody_capture_number_of_reads = 0

def mean(numbers):
    return sum(numbers) / len(numbers) if numbers else 0


for lane, file in file_dict.items():

    print(lane)
    metrics_file = open(file, 'r')

    
    for line in metrics_file.readlines():
        if 'Q30 RNA read' in line:
            
            value1 = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            RNA_Q30_percentages.append(value1)

        if 'Q30 UMI' in line or 'Q30 barcodes' in line:
            value2 = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            RNA_Q30_percentages.append(value2)
    
        if 'Median genes per cell' in line:
            value = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            median_genes_per_cell.append(value)

        if 'Confidently mapped to transcriptome' in line:
            value = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            confidently_mapped_to_transcriptome.append(value)

        if 'Mean antibody reads usable per cell' in line and 'Antibody Capture' in line:
            value = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            mean_antibody_reads_usable_per_cell.append(value)

        if 'Sequencing saturation' in line and 'Gene Expression' in line:
            value = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            gene_expression_sequencing_saturation.append(value)

        if 'Mean used reads per cell' in line and 'VDJ T' in line:
            value = float(line.split(',')[-1].strip().strip('%').strip('"').replace(',',''))
            VDJ_mean_used_reads_per_cell.append(value)

        if 'Number of reads' in line and 'Gene Expression' in line:
            value = float(line.split('"')[-2].strip().strip('%').strip('"').replace(',',''))
            gene_expression_number_of_reads +=  value

        if 'Number of reads' in line and 'VDJ T' in line:
            value = float(line.split('"')[-2].strip().strip('%').strip('"').replace(',',''))
            VDJ_number_of_reads +=  value

        if 'Number of reads' in line and 'Antibody Capture' in line:
            value = float(line.split('"')[-2].strip().strip('%').strip('"').replace(',',''))
            antibody_capture_number_of_reads +=  value

    metrics_file.close()


output_file.write('Range RNA_Q30_percentages \t' + str(min(RNA_Q30_percentages)) + '\t'+ str(max(RNA_Q30_percentages))+ ' Mean ' + str(mean(RNA_Q30_percentages)) + '\n')
output_file.write('Range median_genes_per_cell \t' + str(min(median_genes_per_cell)) + '\t'+ str(max(median_genes_per_cell))+ ' Mean ' + str((mean(median_genes_per_cell))) + '\n')
output_file.write('Range confidently_mapped_to_transcriptome \t' + str(min(confidently_mapped_to_transcriptome)) + '\t'+ str(max(confidently_mapped_to_transcriptome))+ ' Mean '+ str(mean(confidently_mapped_to_transcriptome)) + '\n')
output_file.write('Range mean_antibody_reads_usable_per_cell \t' + str(min(mean_antibody_reads_usable_per_cell)) + '\t'+ str(max(mean_antibody_reads_usable_per_cell))+ ' Mean '+ str(mean(mean_antibody_reads_usable_per_cell)) + '\n') 
output_file.write('Range gene_expression_sequencing_saturation \t' + str(min(gene_expression_sequencing_saturation)) + '\t'+ str(max(gene_expression_sequencing_saturation)) + ' Mean '+ str(mean(gene_expression_sequencing_saturation)) + '\n')
output_file.write('Range VDJ_mean_used_reads_per_cell \t' + str(min(VDJ_mean_used_reads_per_cell)) + '\t'+ str(max(VDJ_mean_used_reads_per_cell))+ ' Mean '+ str(mean(VDJ_mean_used_reads_per_cell)) + '\n')
output_file.write('gene_expression_number_of_reads \t' + str(gene_expression_number_of_reads) + '\n')
output_file.write('VDJ_number_of_reads \t' + str(VDJ_number_of_reads) + '\n')
output_file.write('antibody_capture_number_of_reads \t' + str(antibody_capture_number_of_reads) + '\n')

output_file.close()
