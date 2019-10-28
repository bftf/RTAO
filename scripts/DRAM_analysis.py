import sys
import os
import operator
import matplotlib.pyplot as plt

from enum import IntEnum


'''
======================================================================
THE USER HAS TO CHANGE THESE VARIABLES BASED ON the GPGPU-Sim OUTPUT
======================================================================
'''
'''
filepath = '../build_make/output/teapot_test.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x380

cudaHitsStart = 0xc0000400
cudaHitsSize = 0x1c0

cudaBVHNodeDataStart = 0xc0000600
cudaBVHNodeDataSize = 0x19eab0

cudaInlinedPrimitivesStart = 0xc019f100
cudaInlinedPrimitivesSize = 0x5bafa0 

cudaPrimitiveIndicesStart = 0xc075a100
cudaPrimitiveIndicesSize = 0x16ebe8

start_line = 683
end_line = 13908
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
end_line = 463712
'''

'''
filepath = '../build_make/output/dragon_200k_16_5_oct23.txt'

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

start_line = 688
end_line = 1251209616
'''

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

filepath = '../build_make/output/sanmig_200k_1_5_oct24.txt'

cudaRaysStart = 0xc0000000
cudaRaysSize = 0x61a800

cudaHitsStart = 0xc061a800
cudaHitsSize = 0x30d400

cudaBVHNodeDataStart = 0xc0927c00
cudaBVHNodeDataSize = 0xb06c3d0

cudaInlinedPrimitivesStart = 0xcb994000
cudaInlinedPrimitivesSize = 0x23e88a00 

cudaPrimitiveIndicesStart = 0xef81ca00
cudaPrimitiveIndicesSize = 0x8fa2280

start_line = 687
end_line = 36818413


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

class MemoryLevelEnum(IntEnum):
     LDST_level = 0
     L2_level = 1
     DRAM_level = 2

g_histogram_LDST = [0] * 6
g_histogram_L2 = [0] * 6
g_histogram_DRAM = [0] * 6
g_histogram = [0] * 6

g_reuse_buffer_LDST = {}
g_reuse_buffer_L2 = {}
g_reuse_buffer_DRAM = {}
g_reuse_buffer = {}

def histogram_insert(enum_index, level):
  g_histogram[enum_index] += 1
  if level == MemoryLevelEnum.LDST_level:
    g_histogram_LDST[enum_index] +=1
  elif level == MemoryLevelEnum.L2_level:
    g_histogram_L2[enum_index] +=1
  elif level == MemoryLevelEnum.DRAM_level:
    g_histogram_DRAM[enum_index] +=1

def insert_into_resue_buffer(address, cur_buffer):
  if address in cur_buffer:
    cur_buffer[address] += 1
  else:
    cur_buffer[address] = 1

def reuse_insert(address, enum_index, level):
  insert_into_resue_buffer(address, g_reuse_buffer)
  if level == MemoryLevelEnum.LDST_level:
    insert_into_resue_buffer(address, g_reuse_buffer_LDST)
  elif level == MemoryLevelEnum.L2_level:
    insert_into_resue_buffer(address, g_reuse_buffer_L2)
  elif level == MemoryLevelEnum.DRAM_level:
    insert_into_resue_buffer(address, g_reuse_buffer_DRAM)

def categorize(address, level):
  if address >= cudaRaysStart and address <=cudaRaysStart+cudaRaysSize:
    cur_buffer_enum_index = BufferEnum.cudaRays
  elif address >= cudaHitsStart and address <= cudaHitsStart+cudaHitsSize:
    cur_buffer_enum_index = BufferEnum.cudaHits
  elif address >= cudaBVHNodeDataStart and address <= cudaBVHNodeDataStart+cudaBVHNodeDataSize:
    cur_buffer_enum_index = BufferEnum.cudaBVHNodeData
  elif address >= cudaInlinedPrimitivesStart and address <= cudaInlinedPrimitivesStart+cudaInlinedPrimitivesSize:
    cur_buffer_enum_index = BufferEnum.cudaInlinedPrimitives
  elif address >= cudaPrimitiveIndicesStart and address <= cudaPrimitiveIndicesStart+cudaPrimitiveIndicesSize:
    cur_buffer_enum_index = BufferEnum.cudaPrimitiveIndices
  else:
    cur_buffer_enum_index = BufferEnum.other

  histogram_insert(cur_buffer_enum_index, level)
  reuse_insert(address, cur_buffer_enum_index, level)

def evaluate_histogram():
  print("======== Memory stats ========")
  print("Label: #cudaRays, #cudaHits, #cudaBVHNodeData, #cudaInlinedPrimitives, #cudaPrimitiveIndices, #other")
  print("LDST: ", g_histogram_LDST)
  print("L2: ", g_histogram_L2)
  print("DRAM: ", g_histogram_DRAM)
  print("SUM: ", g_histogram)

def evaluate_reuse():
  print("======== Reuse stats ========")
  for cur_buffer, cur_name in zip([g_reuse_buffer_LDST, g_reuse_buffer_L2, g_reuse_buffer_DRAM, g_reuse_buffer], \
    ["g_reuse_buffer_LDST", "g_reuse_buffer_L2", "g_reuse_buffer_DRAM", "g_reuse_buffer"]):
    squahsed_sum = 0
    total_sum = 0

    for key in cur_buffer:
      squahsed_sum += 1
      total_sum += cur_buffer[key]
    
    print("reuse stats for: ", cur_name)
    print("ideal: ", squahsed_sum)
    print("real: ", total_sum)
    print("ratio: ", squahsed_sum / total_sum)
    print("inv - ratio: ", total_sum / squahsed_sum)


def main():

  if not os.path.isfile(filepath):
     print("File path {} does not exist. Exiting...".format(filepath))
     sys.exit()

  memory_accesses = {}

  with open(filepath) as fp:
    cur_line_count = 0
    for line in fp:
      cur_line_count += 1
      if cur_line_count >= end_line:
        break
      elif cur_line_count < start_line:
        continue
      else:
        line = line.replace(':', ',')
        cur_list = line.split(",")
        
        if len(cur_list) != 5:
          print("Problem parsing this memory access (skipping this line):")
          print(cur_list)
          continue

        level_label = cur_list[0]
        sid = int(cur_list[1])
        wip = int(cur_list[2])
        address = int(cur_list[3], 16) # interpret as hex


        if level_label == "LDST":
          cur_label =  MemoryLevelEnum.LDST_level
        elif level_label == "L2":
          cur_label =  MemoryLevelEnum.L2_level
        elif level_label == "DRAM":
          cur_label =  MemoryLevelEnum.DRAM_level
        else:
          print("Problem categorizing this memory access (skipping this line):")
          print(cur_list)
          continue
        
        categorize(address, cur_label)



    evaluate_histogram()
    evaluate_reuse()
    

if __name__ == '__main__':
  main()