import sys
import os
import operator

from enum import IntEnum


'''
======================================================================
THE USER HAS TO CHANGE THESE VARIABLES BASED ON the GPGPU-Sim OUTPUT
======================================================================
'''

# filepath = '../build_make/output/sponza_200k_rayFile_nov7_bypassL2.txt'
# filepath = '../build_make/output/debugging.txt'
# filepath = '../build_make/output/sponza_200k_rayFile_nov12_origin_nodf.txt'
filepath =  '../build_make/output/sponza_200k_rayFile_nov13_random_nodf.txt'


'''
======================================================================
'''

cudaRaysStart = 0x0
cudaRaysSize = 0x0 
cudaHitsStart = 0x0
cudaHitsSize = 0x0
cudaBVHNodeDataStart = 0x0
cudaBVHNodeDataSize = 0x0
cudaInlinedPrimitivesStart = 0x0
cudaInlinedPrimitivesSize = 0x0 
cudaPrimitiveIndicesStart = 0x0
cudaPrimitiveIndicesSize = 0x0

class BufferEnum(IntEnum):
     cudaRays = 0
     cudaHits = 1
     cudaBVHNodeData = 2
     cudaInlinedPrimitives = 3
     cudaPrimitiveIndices = 4
     other = 5
     instructions = 6
     num_elements = 7

class MemoryLevelEnum(IntEnum):
     LDST_level = 0
     L2_level = 1
     DRAM_level = 2

     L2_hit = 3
     L2_miss = 4
     L2_rf = 5

     L1C_hit = 6
     L1C_miss = 7
     L1C_rf = 8

     L1T_hit = 9
     L1T_miss = 10
     L1T_rf = 11

     L1D_hit = 12
     L1D_miss = 13
     L1D_rf = 14

     L1I_hit = 15
     L1I_miss = 16
     L1I_rf = 17

g_histogram = [0] * BufferEnum.num_elements

g_histogram_LDST = [0] * BufferEnum.num_elements
g_histogram_L2 = [0] * BufferEnum.num_elements
g_histogram_DRAM = [0] * BufferEnum.num_elements

g_histogram_L2_hit = [0] * BufferEnum.num_elements
g_histogram_L2_miss = [0] * BufferEnum.num_elements
g_histogram_L2_rf = [0] * BufferEnum.num_elements

g_histogram_L1C_hit = [0] * BufferEnum.num_elements
g_histogram_L1C_miss = [0] * BufferEnum.num_elements
g_histogram_L1C_rf = [0] * BufferEnum.num_elements

g_histogram_L1T_hit = [0] * BufferEnum.num_elements
g_histogram_L1T_miss = [0] * BufferEnum.num_elements
g_histogram_L1T_rf = [0] * BufferEnum.num_elements

g_histogram_L1D_hit = [0] * BufferEnum.num_elements
g_histogram_L1D_miss = [0] * BufferEnum.num_elements
g_histogram_L1D_rf = [0] * BufferEnum.num_elements

g_histogram_L1I_hit = [0] * BufferEnum.num_elements
g_histogram_L1I_miss = [0] * BufferEnum.num_elements
g_histogram_L1I_rf = [0] * BufferEnum.num_elements

g_reuse_buffer_LDST = {}
g_reuse_buffer_L2 = {}
g_reuse_buffer_DRAM = {}
g_reuse_buffer = {}

def histogram_insert(enum_index, level, cur_buffer):
  g_histogram[enum_index] += 1
  cur_buffer[enum_index] += 1

def insert_into_resue_buffer(address, cur_buffer):
  if address in cur_buffer:
    cur_buffer[address] += 1
  else:
    cur_buffer[address] = 1

def categorize(address):
  if address >= cudaRaysStart and address < cudaHitsStart:
    cur_buffer_enum_index = BufferEnum.cudaRays
  elif address >= cudaHitsStart and address < cudaBVHNodeDataStart:
    cur_buffer_enum_index = BufferEnum.cudaHits
  elif address >= cudaBVHNodeDataStart and address < cudaInlinedPrimitivesStart:
    cur_buffer_enum_index = BufferEnum.cudaBVHNodeData
  elif address >= cudaInlinedPrimitivesStart and address < cudaPrimitiveIndicesStart:
    cur_buffer_enum_index = BufferEnum.cudaInlinedPrimitives
  elif address >= cudaPrimitiveIndicesStart and address < cudaPrimitiveIndicesStart+cudaPrimitiveIndicesSize:
    cur_buffer_enum_index = BufferEnum.cudaPrimitiveIndices
  elif address >= 0xf0000000:
    cur_buffer_enum_index = BufferEnum.instructions
  else:
    cur_buffer_enum_index = BufferEnum.other


  return cur_buffer_enum_index
  

def evaluate_histogram():
  print("======== Memory stats ========")
  print("Label: #cudaRays, #cudaHits, #cudaBVHNodeData, #cudaInlinedPrimitives, #cudaPrimitiveIndices, #other, #instructions")
  print("LDST: ", g_histogram_LDST)
  print("L2: ", g_histogram_L2)
  print("DRAM: ", g_histogram_DRAM)

  
  print("L2_hit: ", g_histogram_L2_hit)
  print("L2_miss: ", g_histogram_L2_miss)
  print("L2_rf: ", g_histogram_L2_rf)

  print("L1D_hit: ", g_histogram_L1D_hit)
  print("L1D_miss: ", g_histogram_L1D_miss)
  print("L1D_rf: ", g_histogram_L1D_rf)

  print("L1C_hit: ", g_histogram_L1C_hit)
  print("L1C_miss: ", g_histogram_L1C_miss)
  print("L1C_rf: ", g_histogram_L1C_rf)

  print("L1T_hit: ", g_histogram_L1T_hit)
  print("L1T_miss: ", g_histogram_L1T_miss)
  print("L1T_rf: ", g_histogram_L1T_rf)

  print("L1I_hit: ", g_histogram_L1I_hit)
  print("L1I_miss: ", g_histogram_L1I_miss)
  print("L1I_rf: ", g_histogram_L1I_rf)
  

  print("SUM: ", g_histogram)

def evaluate_reuse_in_buffer(cur_buffer, name):
  
  total_number_of_addresses = len(cur_buffer)
  print('Total number of accesses in ' + name + ': ' + str(total_number_of_addresses))

  if total_number_of_addresses == 0:
    print("total_number_of_addresses == 0, skipping buffer: " + name)
    return 

  i = 0
  total_accesses = 0

  top_1_sum = 0
  top_5_sum = 0
  top_10_sum = 0
  top_20_sum = 0
  top_30_sum = 0
  top_40_sum = 0

  top_1_cutoff = -1
  top_5_cutoff = -1
  top_10_cutoff = -1
  top_20_cutoff = -1
  top_30_cutoff = -1
  top_40_cutoff = -1

  top_1_categories = [0] * BufferEnum.num_elements
  top_5_categories = [0] * BufferEnum.num_elements
  top_10_categories = [0] * BufferEnum.num_elements
  top_20_categories = [0] * BufferEnum.num_elements
  top_30_categories = [0] * BufferEnum.num_elements
  top_40_categories = [0] * BufferEnum.num_elements
  
  cur_buffer_sorted = sorted(cur_buffer.items(), key=operator.itemgetter(1), reverse=True)

  for address, num_access in cur_buffer_sorted:
    if i < 0.01 * total_number_of_addresses:
      top_1_sum += num_access
      top_1_categories[categorize(address)] += 1
      top_1_cutoff = num_access
    if i < 0.05 * total_number_of_addresses:
      top_5_sum += num_access
      top_5_categories[categorize(address)] += 1
      top_5_cutoff = num_access
    if i < 0.1 * total_number_of_addresses:
      top_10_sum += num_access
      top_10_categories[categorize(address)] += 1
      top_10_cutoff = num_access
    if i < 0.2 * total_number_of_addresses:
      top_20_sum += num_access
      top_20_categories[categorize(address)] += 1
      top_20_cutoff = num_access
    if i < 0.3 * total_number_of_addresses:
      top_30_sum += num_access
      top_30_categories[categorize(address)] += 1
      top_30_cutoff = num_access
    if i < 0.4 * total_number_of_addresses:
      top_40_sum += num_access
      top_40_categories[categorize(address)] += 1
      top_40_cutoff = num_access

    total_accesses += num_access
    i += 1

  print("reuse stats for: " + name)
  print("top 1%: " + str(top_1_sum / total_accesses))
  print("top 5%: " + str(top_5_sum / total_accesses))
  print("top 10%: " + str(top_10_sum / total_accesses))
  print("top 20%: " + str(top_20_sum / total_accesses))
  print("top 30%: " + str(top_30_sum / total_accesses))
  print("top 40%: " + str(top_40_sum / total_accesses))

  print("top 1% categories: " + str(top_1_categories))
  print("top 5% categories: " + str(top_5_categories))
  print("top 10% categories: " + str(top_10_categories))
  print("top 20% categories: " + str(top_20_categories))
  print("top 30% categories: " + str(top_30_categories))
  print("top 40% categories: " + str(top_40_categories))

  print("top 1% cutoff: " + str(top_1_cutoff))
  print("top 5% cutoff: " + str(top_5_cutoff))
  print("top 10% cutoff: " + str(top_10_cutoff))
  print("top 20% cutoff: " + str(top_20_cutoff))
  print("top 30% cutoff: " + str(top_30_cutoff))
  print("top 40% cutoff: " + str(top_40_cutoff))


def evaluate_reuse():
  evaluate_reuse_in_buffer(g_reuse_buffer_LDST, 'LDST')
  evaluate_reuse_in_buffer(g_reuse_buffer_L2, 'L2')
  evaluate_reuse_in_buffer(g_reuse_buffer_DRAM, 'DRAM')

def evaluate_hit_rate(name, hit_buffer, miss_buffer):
  miss_rates = [0] * (BufferEnum.num_elements + 1)
  total_miss = 0
  total_access = 0
  
  for index, (cur_hit_count, cur_miss_count) in enumerate(zip(hit_buffer, miss_buffer)):
    total_miss += cur_miss_count
    total_access += cur_miss_count + cur_hit_count
    miss_rates[index] = cur_miss_count / (cur_miss_count + cur_hit_count) if cur_miss_count else 0
    
  miss_rates[BufferEnum.num_elements] = total_miss / total_access if total_access else 0
  print("Miss rate " + name +": " + str(miss_rates))

def evaluate_hit_rates():
  print("Label: #cudaRays, #cudaHits, #cudaBVHNodeData, #cudaInlinedPrimitives, #cudaPrimitiveIndices, #other, #instructions, total")
  evaluate_hit_rate("L1D", g_histogram_L1D_hit, g_histogram_L1D_miss)
  evaluate_hit_rate("L1C", g_histogram_L1C_hit, g_histogram_L1C_miss)
  evaluate_hit_rate("L2", g_histogram_L2_hit, g_histogram_L2_miss)

l2_addresses = {}
l2_line_size = 128
cold_miss = 0
l2_hits = 0

def insert_l2_address(addr):
  global cold_miss
  global l2_hits 
  global l2_line_size

  block_address = addr & ~(l2_line_size-1)

  if block_address not in l2_addresses:
    l2_addresses[block_address] = 1
    cold_miss += 1
  else:
    l2_addresses[block_address] += 1
    l2_hits += 1

dram_addresses = {}

def insert_dram_address(addr):
  global l2_line_size

  block_address = addr & ~(l2_line_size-1)

  if block_address not in dram_addresses:
    dram_addresses[block_address] = 1
  else:
    dram_addresses[block_address] += 1

def compute_l2_cold_miss_rate():
  print("l2 cold misses: ", cold_miss)
  print("l2 max hits: ", l2_hits)
  print("l2 upper limit miss rate: ", cold_miss/(l2_hits + cold_miss))
  print("actual DRAM accesses: ", sum(g_histogram_DRAM))

  sorted_l2_addresses = sorted(l2_addresses.items(), key=operator.itemgetter(1), reverse=True)
  
  for index, (address, hit_count) in enumerate(sorted_l2_addresses):
    # do lookup in DRAM histogram
    dram_count = 0
    if address in dram_addresses:
      dram_count = dram_addresses[address]

    print(str(hex(address)) + " was hit: " + str(hit_count) + ", label: " + str(categorize(address)) + ", dram count: " + str(dram_count))
    if index > 1000:
      break

def main():

  global cudaRaysStart
  global cudaRaysSize 
  global cudaHitsStart
  global cudaHitsSize
  global cudaBVHNodeDataStart
  global cudaBVHNodeDataSize
  global cudaInlinedPrimitivesStart
  global cudaInlinedPrimitivesSize 
  global cudaPrimitiveIndicesStart
  global cudaPrimitiveIndicesSize

  if not os.path.isfile(filepath):
     print("File path {} does not exist. Exiting...".format(filepath))
     sys.exit()

  memory_accesses = {}
  first_mem_access = True

  with open(filepath) as fp:
    cur_line_count = 0
    for line in fp:
      cur_line_count += 1

      if "Destroy streams for kernel" in line:
        print("Exit keywords detected")
        break
      elif "cudaRays" in line and len(line.split(" ")) == 3:
        cur_list = line.split(" ")
        cudaRaysStart = int(cur_list[1], 16)
        cudaRaysSize = int(cur_list[2], 16)
      elif "cudaHits" in line and len(line.split(" ")) == 3:
        cur_list = line.split(" ")
        cudaHitsStart = int(cur_list[1], 16)
        cudaHitsSize = int(cur_list[2], 16)
      elif "cudaBVHNodeData" in line and len(line.split(" ")) == 3:
        cur_list = line.split(" ")
        cudaBVHNodeDataStart = int(cur_list[1], 16)
        cudaBVHNodeDataSize = int(cur_list[2], 16)
      elif "cudaInlinedPrimitives" in line and len(line.split(" ")) == 3:
        cur_list = line.split(" ")
        cudaInlinedPrimitivesStart = int(cur_list[1], 16)
        cudaInlinedPrimitivesSize = int(cur_list[2], 16)
      elif "cudaPrimitiveIndices" in line and len(line.split(" ")) == 3:
        cur_list = line.split(" ")
        cudaPrimitiveIndicesStart = int(cur_list[1], 16)
        cudaPrimitiveIndicesSize = int(cur_list[2], 16)
      else:
        line = line.replace(':', ',')
        cur_list = line.split(" ")
        
        if len(cur_list) != 6:
          # print("Problem parsing this memory access (skipping this line):")
          # print(cur_list)
          continue
        try:
          level_label = cur_list[0]
          sid = int(cur_list[1], 16) # interpret as hex
          wip = int(cur_list[2], 16)
          address = int(cur_list[3], 16)
          read_or_write = int(cur_list[4]) 
          size = int(cur_list[5], 16)
        except:
          print(cur_list)
          continue

        if (first_mem_access == True):
          first_mem_access = False
          assert(cudaRaysStart > 0 and cudaRaysSize > 0 and cudaBVHNodeDataStart > 0 and cudaBVHNodeDataSize > 0)
          assert(cudaInlinedPrimitivesStart > 0 and cudaInlinedPrimitivesSize > 0 and cudaPrimitiveIndicesStart > 0 and cudaPrimitiveIndicesSize > 0)

          print('cudaRaysStart: ' + str(hex(cudaRaysStart)))
          print('cudaRaysSize: ' + str(hex(cudaRaysSize)))
          print('cudaHitsStart: ' + str(hex(cudaHitsStart)))
          print('cudaHitsSize: ' + str(hex(cudaHitsSize)))
          print('cudaBVHNodeDataStart: ' + str(hex(cudaBVHNodeDataStart)))
          print('cudaBVHNodeDataSize: ' + str(hex(cudaBVHNodeDataSize)))
          print('cudaInlinedPrimitivesStart: ' + str(hex(cudaInlinedPrimitivesStart)))
          print('cudaInlinedPrimitivesSize: ' + str(hex(cudaInlinedPrimitivesSize)))
          print('cudaPrimitiveIndicesStart: ' + str(hex(cudaPrimitiveIndicesStart))) 
          print('cudaPrimitiveIndicesSize: ' + str(hex(cudaPrimitiveIndicesSize)))
          print('Start Parsing started')

        cur_reuse = -1

        if level_label == "LDST":
          # this is no longer reached
          assert(False)
        elif level_label == "L2":
          # this is no longer reached
          assert(False)
        elif level_label == "DRAM":
          cur_label =  MemoryLevelEnum.DRAM_level
          cur_hist = g_histogram_DRAM
          cur_reuse = g_reuse_buffer_DRAM
        elif level_label == "L2HIT":
          cur_label =  MemoryLevelEnum.L2_hit
          cur_hist = g_histogram_L2_hit
          cur_reuse = g_reuse_buffer_L2
        elif level_label == "L2MISS":
          cur_label =  MemoryLevelEnum.L2_miss
          cur_hist = g_histogram_L2_miss
          cur_reuse = g_reuse_buffer_L2
        elif level_label == "L2RF":
          cur_label =  MemoryLevelEnum.L2_rf
          cur_hist = g_histogram_L2_rf
        elif level_label == "L1IHIT":
          cur_label =  MemoryLevelEnum.L1I_hit
          cur_hist = g_histogram_L1I_hit
        elif level_label == "L1IMISS":
          cur_label =  MemoryLevelEnum.L1I_miss
          cur_hist = g_histogram_L1I_miss
        elif level_label == "L1IRF":
          cur_label =  MemoryLevelEnum.L1I_rf
          cur_hist = g_histogram_L1I_rf
        elif level_label == "L1DHIT":
          cur_label =  MemoryLevelEnum.L1D_hit
          cur_hist = g_histogram_L1D_hit
        elif level_label == "L1DMISS":
          cur_label =  MemoryLevelEnum.L1D_miss
          cur_hist = g_histogram_L1D_miss
        elif level_label == "L1DRF":
          cur_label =  MemoryLevelEnum.L1D_rf
          cur_hist = g_histogram_L1D_rf
        elif level_label == "L1CHIT":
          cur_label =  MemoryLevelEnum.L1C_hit
          cur_hist = g_histogram_L1C_hit
        elif level_label == "L1CMISS":
          cur_label =  MemoryLevelEnum.L1C_miss
          cur_hist = g_histogram_L1C_miss
        elif level_label == "L1CRF":
          cur_label =  MemoryLevelEnum.L1C_rf
          cur_hist = g_histogram_L1C_rf
        elif level_label == "L1THIT":
          cur_label =  MemoryLevelEnum.L1T_hit
          cur_hist = g_histogram_L1T_hit
        elif level_label == "L1TMISS":
          cur_label =  MemoryLevelEnum.L1T_miss
          cur_hist = g_histogram_L1T_miss
        elif level_label == "L1TRF":
          cur_label =  MemoryLevelEnum.L1T_rf
          cur_hist = g_histogram_L1T_rf
        else:
          print("Problem categorizing this memory access (skipping this line):")
          print(cur_list)
          continue

        #error checking - make sure all instructions have address higher than 0xf0000000
        if level_label == 'L1IHIT' or level_label == 'L1IMISS':
          assert(address >= 0xf0000000)
        elif level_label == 'L1DHIT' or level_label == 'L1DMISS' \
          or level_label == 'L1THIT' or level_label == 'L1TMISS' \
          or level_label == 'L1CHIT' or level_label == 'L1CMISS':
            assert(address < 0xf0000000)

        # primary classification
        cur_buffer_enum_index = categorize(address)
        histogram_insert(cur_buffer_enum_index, cur_label, cur_hist)

        # secondary classification - if an access was for any cache on the L1 level, 
        # it should be added to the LDST as well!
        if cur_label >= 6:
          cur_label = MemoryLevelEnum.LDST_level
          cur_hist = g_histogram_LDST
          histogram_insert(cur_buffer_enum_index, cur_label, cur_hist)
          cur_reuse = g_reuse_buffer_LDST
        
        # filter instructions out when evaluating reuse
        if cur_reuse is not -1 and cur_buffer_enum_index is not BufferEnum.instructions:
          insert_into_resue_buffer(address, cur_reuse)

        # compute l2 cold misses
        if level_label == 'L2HIT' or level_label == 'L2MISS':
          insert_l2_address(address)

        # generate DRAM histogram
        if level_label == 'DRAM':
          insert_dram_address(address)

    # print stats
    evaluate_histogram()
    evaluate_reuse()
    evaluate_hit_rates()
    compute_l2_cold_miss_rate()

if __name__ == '__main__':
  main()
