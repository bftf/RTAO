import sys
import os
import operator
import matplotlib.pyplot as plt

def main():
   
  #filepath = '../build_make/output/cwbvh_sponza_1_1_10k_DRAMExperiment.txt'
  filepath = '../build_make/output/cwbvh_dragon_1_4_50k_DRAMExperiment.txt'
  
  start_line = 678
  end_line = 220070

  if not os.path.isfile(filepath):
     print("File path {} does not exist. Exiting...".format(filepath))
     sys.exit()

  memory_accesses = {}
  discarded_accesses_count = 0

  with open(filepath) as fp:
    cur_line_count = 0
    for line in fp:
      cur_line_count += 1
      if cur_line_count >= end_line:
        break
      elif cur_line_count < start_line:
        continue
      else:
        cur_list = line.split(",")
        
        if len(cur_list) != 5:
          print("Problem parsing this memory access (skipping this line):")
          print(cur_list)
          continue

        sid = int(cur_list[0])
        wip = int(cur_list[1])
        cur_dram_access = cur_list[2]


        if sid < 1000:
          if cur_dram_access in memory_accesses:
            memory_accesses[cur_dram_access] += 1
          else:
            memory_accesses[cur_dram_access] = 1
        else:
          discarded_accesses_count += 1

    print(discarded_accesses_count)

    sorted_memory_accesses = sorted(memory_accesses.items(), key=operator.itemgetter(1), reverse=True)

    for i in range(0, 100):
      print(hex(int(sorted_memory_accesses[i][0], 16)), sorted_memory_accesses[i])

    squahsed_sum = 0
    total_sum = 0

    for cur_access in sorted_memory_accesses:
      squahsed_sum += 1
      total_sum += cur_access[1]
    print("ideal: ", squahsed_sum)
    print("real: ", total_sum)
    print("ratio: ", squahsed_sum / total_sum)
    print("inv - ratio: ", total_sum / squahsed_sum)

    

if __name__ == '__main__':
  main()