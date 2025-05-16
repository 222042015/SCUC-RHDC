# iterate instance name in the list
threads=4
nthreads=4

timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo "==========Running script on instance: $timestamp, threads: $threads, nThreads: $nthreads ==========" >> res_an_threads_${threads}_nthreads_${nthreads}.txt
for instance in $(cat instances.txt); do
  julia --threads=$threads test_cb.jl -dataset $instance -nInt 6 -nCont 6 -stepsize 6 -threads $nthreads | tee ./logs/rh/${instance}_threads_${threads}_nthreads_${nthreads}.log;
done

timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo "==========Running script on instance: $timestamp, threads: $threads, nThreads: $nthreads ==========" >> res_td_threads_${threads}_nthreads_${nthreads}.txt
for instance in $(cat instances.txt); do
  julia --threads=$threads test_baseline.jl -dataset $instance -threads $nthreads| tee ./logs/td/${instance}_threads_${threads}_nthreads_${nthreads}.log;
done


timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo "==========Running script on instance: $timestamp, nThreads: $nthreads ==========" >> res_uc_threads_${threads}_nthreads_${nthreads}.txt
for instance in $(cat instances.txt); do
  julia --threads=$threads test_uc.jl -dataset $instance -threads $nthreads| tee ./logs/uc/${instance}_threads_${threads}_nthreads_${nthreads}.log;
done