while read sample;
	do qsub job/${sample}.sh;
done < pending.vtr
