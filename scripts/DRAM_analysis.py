import sys
import os
import operator

from enum import IntEnum


'''
======================================================================
THE USER HAS TO CHANGE THESE VARIABLES BASED ON the GPGPU-Sim OUTPUT
======================================================================
'''
'''
filepath = '../build_make/output/teapot_test.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0xc80 # 100 rays

cudaHitsStart = 0xc0000d00
cudaHitsSize = 0x640

cudaBVHNodeDataStart = 0xc0001400
cudaBVHNodeDataSize = 0x177a0

cudaInlinedPrimitivesStart = 0xc0018c00
cudaInlinedPrimitivesSize = 0x6a620 

cudaPrimitiveIndicesStart = 0xc0083300
cudaPrimitiveIndicesSize = 0x1a988

start_line = 647
end_line = 100000
'''

'''
filepath = '../build_make/output/dragon_50k_10_23.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x27100

cudaHitsStart = 0xc0027100
cudaHitsSize = 0x13880

cudaBVHNodeDataStart = 0xc003aa00
cudaBVHNodeDataSize = 0x19eab0

cudaInlinedPrimitivesStart = 0xc01d9500
cudaInlinedPrimitivesSize = 0x5bafa0 

cudaPrimitiveIndicesStart = 0xc0794500
cudaPrimitiveIndicesSize = 0x16ebe8

start_line = 688
end_line = 5000
'''

filepath = '../build_make/output/dragon_200k_rayFile_nov3_bypassL2.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x61a800

cudaHitsStart = 0xc061a800
cudaHitsSize = 0x30d400

cudaBVHNodeDataStart = 0xc0927c00
cudaBVHNodeDataSize = 0x19eab0

cudaInlinedPrimitivesStart = 0xc0ac6700
cudaInlinedPrimitivesSize = 0x5bafa0 

cudaPrimitiveIndicesStart = 0xc1081700
cudaPrimitiveIndicesSize = 0x16ebe8

start_line = 647
end_line = 60000000

'''
filepath = '../build_make/output/sponza_200k_16_5_oct23.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x61a800

cudaHitsStart = 0xc061a800
cudaHitsSize = 0x30d400

cudaBVHNodeDataStart = 0xc0927c00
cudaBVHNodeDataSize = 0x4ac9f0

cudaInlinedPrimitivesStart = 0xc0dd4600
cudaInlinedPrimitivesSize = 0x11c5590 

cudaPrimitiveIndicesStart = 0xc1f99c00
cudaPrimitiveIndicesSize = 0x471564

start_line = 1880
end_line = 12994282
'''

'''
filepath = '../build_make/output/teapot_200k_16_5_oct24.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x61a800

cudaHitsStart = 0xc061a800
cudaHitsSize = 0x30d400

cudaBVHNodeDataStart = 0xc0927c00
cudaBVHNodeDataSize = 0x177a0

cudaInlinedPrimitivesStart = 0xc093f400
cudaInlinedPrimitivesSize = 0x6a620 

cudaPrimitiveIndicesStart = 0xc09a9b00
cudaPrimitiveIndicesSize = 0x1a988

start_line = 688
end_line = 10058076
'''

'''
filepath = '../build_make/output/sanmig_200k_rayFile_L2_10_30.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x61a800

cudaHitsStart = 0xc061a800
cudaHitsSize = 0x30d400

cudaBVHNodeDataStart = 0xc0927c00
cudaBVHNodeDataSize = 0x10452bc0

cudaInlinedPrimitivesStart = 0xd0d7a800
cudaInlinedPrimitivesSize = 0x332a2790 

cudaPrimitiveIndicesStart = 0x10401d000
cudaPrimitiveIndicesSize = 0xcca89e4

start_line = 686
end_line = 3393038
'''


'''
======================================================================
'''

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
  if address >= cudaRaysStart and address < cudaRaysStart+cudaRaysSize:
    cur_buffer_enum_index = BufferEnum.cudaRays
  elif address >= cudaHitsStart and address < cudaHitsStart+cudaHitsSize:
    cur_buffer_enum_index = BufferEnum.cudaHits
  elif address >= cudaBVHNodeDataStart and address < cudaBVHNodeDataStart+cudaBVHNodeDataSize:
    cur_buffer_enum_index = BufferEnum.cudaBVHNodeData
  elif address >= cudaInlinedPrimitivesStart and address < cudaInlinedPrimitivesStart+cudaInlinedPrimitivesSize:
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

def main():

  if not os.path.isfile(filepath):
     print("File path {} does not exist. Exiting...".format(filepath))
     sys.exit()

  memory_accesses = {}

  with open(filepath) as fp:
    cur_line_count = 0
    for line in fp:
      cur_line_count += 1

      if "Destroy streams for kernel" in line:
        print("Exit keywords detected")
        break
      if cur_line_count >= end_line:
        break
      elif cur_line_count < start_line:
        continue
      else:
        line = line.replace(':', ',')
        cur_list = line.split(",")
        
        if len(cur_list) != 5:
          #print("Problem parsing this memory access (skipping this line):")
          #print(cur_list)
          continue

        level_label = cur_list[0]
        sid = int(cur_list[1])
        wip = int(cur_list[2])
        address = int(cur_list[3], 16) # interpret as hex

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


    # print stats
    evaluate_histogram()
    evaluate_reuse()
    

if __name__ == '__main__':
  main()
