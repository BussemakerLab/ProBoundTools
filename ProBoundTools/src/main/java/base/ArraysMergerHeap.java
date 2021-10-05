package base;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.LinkedList;
import java.util.Properties;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

/**
 * A heap structure to merge CountObject files. Each node in the heap is mapped to one file.
 */
public class ArraysMergerHeap
{
	public static final int BUFFER_SIZE_IN  = 1024*1024 * 200; //200M in total
	public static final int BUFFER_SIZE_OUT = 1024*1024 * 30;  //30M in total
	private int BUFFER_SIZE_PER_FILE =  1024*1024 ;
	private MinHeapNode[] arrays;
	private int len;
	
	public ArraysMergerHeap(String tempFolder, String filePrefix, int size)
	{
		this.len= size;
		this.arrays = new MinHeapNode[len];

		int MAX_THREAD_NUM =  Util.getMaxThreadNumber();
		BUFFER_SIZE_PER_FILE =  BUFFER_SIZE_IN / size; 
		DebugLog.log("Buffer size per file:"+ BUFFER_SIZE_PER_FILE);
		int total_quque_size = 1024*1024*MAX_THREAD_NUM*2;
		DebugLog.log("total_quque_size : "+ total_quque_size);
		int QUEUE_SIZE = (int) (total_quque_size / size);
		
		for(int i=0;i<size;i++)
		{
			this.arrays[i]=new MinHeapNode( ArraysMerger.getFileName(tempFolder, filePrefix, i),
					BUFFER_SIZE_PER_FILE, QUEUE_SIZE);
			//Util.printMemeoryUsage();
		}
		
	}
	
	public static class MinHeapNode
	{
		CountObject element;
		boolean end=false;
		String path;
		long currentPosition = 0;
		LinkedList<CountObject> queue;
		final int OBJECT_SIZE=16;
		int QUEUE_SIZE=100000;
		int BUFFER_SIZE_PER_FILE = 1024*1024;
		
		public MinHeapNode(String path, int bufferSize, int queueSize)
		{
			this.path=path;
			this.BUFFER_SIZE_PER_FILE = bufferSize;
			this.QUEUE_SIZE = queueSize ;

			DebugLog.log("QUEUE_SIZE:"+this.QUEUE_SIZE);
			DebugLog.log("Data file path:"+this.path);
			reset();
		}

		public MinHeapNode(String path)
		{
			this.path=path;
			DebugLog.log("Data file path:"+this.path);
			reset();
		}
		
		public void reset()
		{
			try
			{
				currentPosition = 0;
				end=false;
			} catch (Exception e)
			{
				DebugLog.log(e);
			}
		}
		
		public CountObject peek()
		{
			if(end)
				return null;
			if(element==null)
			{
				read();
			}
			return element;
		}
		
		public void pop()
		{
			if(!end)
				read();
		}
		
		private void readInElements()
		{
			//DebugLog.log("Reading in new elements from "+path);
			queue=new LinkedList<CountObject>();
			try
			{
				DataInputStream in=new DataInputStream(new BufferedInputStream(new FileInputStream(path),
						BUFFER_SIZE_PER_FILE));
				long skippedBytes = in.skip(currentPosition);
				//DebugLog.log("Skipping to position:"+currentPosition);
				if(skippedBytes<currentPosition)
				{
					return;
				}

				//Util.printMemeoryUsage();
				for(int i=0;i<QUEUE_SIZE;i++)
				{
					element=CountObject.steamInSequence(in);
					if(element==null)
					{
						break;
					}
					queue.addLast(element);
				}
				currentPosition += 16*queue.size();
				//DebugLog.log("Number of element read:"+queue.size());
				//DebugLog.log("Current position:"+currentPosition);
				in.close();
				
				//Util.printMemeoryUsage();
			} catch (Exception e)
			{
				DebugLog.log(e);
				throw new RuntimeException(e);
			}
		}
		
		private void read()
		{
			if(end)
			{
				element=null;
				return;
			}
			try
			{
				if(this.queue!=null && this.queue.size()>0)
				{
					element=this.queue.removeFirst();
				}
				else
				{
					//DebugLog.log("Renewing buffer "+ this.path);
					readInElements(); //renew the buffer
					if(this.queue.size()==0)
					{
						end=true;
						element = null;
					}
					else
					{
						element=this.queue.removeFirst();
					}
				}
			} 
			catch (Exception e)
			{
				end=true;
				// TODO Auto-generated catch block
				DebugLog.log(e);
				System.err.println("Error file:"+path);
				throw new RuntimeException(e);
			} 
		}
	}
	
	private boolean isLess(CountObject o1, CountObject o2)
	{
		if(o1==null) return false;
		if(o1.compareToByKey(o2)<0)return true;
		return false;
	}
	
	public void minHeapify(int i)
	{
		int l=2*i;
		int r=l+1;
		int smallest = i;
		
		CountObject i_ele = arrays[i].peek();
		if(l<len )
		{
			CountObject l_ele = arrays[l].peek();
			if(l_ele==null)
				smallest = i;
			else if(i_ele==null)
				smallest = l;
			else if(isLess(l_ele,i_ele))
			{
				smallest=l;
			}
		}
		CountObject s_ele = arrays[smallest].peek();
		if(r<len )
		{
			CountObject r_ele = arrays[r].peek();
			if(r_ele==null)
			{}
			else if(s_ele==null)
				smallest = r;
			else if(isLess(r_ele,s_ele))
			{
				smallest=r;
			}
		}
		
		if(smallest!=i && arrays[smallest].peek()!=null)
		{
			exchange(smallest,i);
			minHeapify(smallest);
		}
	}
	
	public void buildMinHeap()
	{
		for(int i=(len-1)/2;i>=0;i--)
		{
			//DebugLog.log(Arrays.toString(array));
			minHeapify(i);
		}
	}
	
	/**
	 * Finds the smallest element in all the available file, pops it, loads new element, 
	 * 		and adjusts the heap accordingly.
	 * 
	 * Repeats the above step until all the elements are null.
	 * 
	 * @param outputText
	 * @param filePath
	 * @return
	 * @throws Exception
	 */
	public Properties output(boolean outputText, String filePath) throws Exception
	{
		//DebugLog.log(Arrays.toString(array));
		buildMinHeap();
		
		CountObject obj=null;
		long size=0;
		
		DebugLog.log("Outputing to "+filePath);
		DataOutputStream out = new DataOutputStream(
				new BufferedOutputStream(new FileOutputStream(filePath), BUFFER_SIZE_OUT));
		
		long totalCounts=0;
		long lowestCount=Long.MAX_VALUE;
		long highestCount =Long.MIN_VALUE;
		long MSG_ROWS = 2000000;
		while(true)
		{
			CountObject temp=arrays[0].peek();
			if(temp==null)
			{
				break;
			}
			//DebugLog.log("Size:"+size);
			
			if(obj==null)
			{
				obj=temp;
			}
			else
			{
				if(obj.getKey().compareTo(temp.getKey())==0)
				{
					obj.setCount(obj.getCount() + temp.getCount());
				}
				else
				{
					//DebugLog.log(obj); //output
					if(outputText)
					{
						CountObject.steamOutSequenceTxt(out,obj);
					}
					else
					{
						CountObject.steamOutSequence(out,obj);
					}
					totalCounts+=obj.getCount();
					if (lowestCount > obj.getCount())
					{
						lowestCount = obj.getCount();
					}
					if (highestCount < obj.getCount())
					{
						highestCount = obj.getCount();
					}
					
					//
					size++;
					obj = temp;
					//
					if(size == MSG_ROWS)
					{
						DebugLog.log("Output Size:"+ size);
						MSG_ROWS += 2000000;
					}
				}
			}
			
			//DebugLog.log(arrays[0].peek());
			arrays[0].pop();
			minHeapify(0);
		}

		//DebugLog.log(obj);
		if(obj!=null)
		{
			if(outputText)
			{
				CountObject.steamOutSequenceTxt(out,obj);
			}
			else
			{
				CountObject.steamOutSequence(out,obj);
			}
			totalCounts+=obj.getCount();
			if (lowestCount > obj.getCount())
			{
				lowestCount = obj.getCount();
			}
			if (highestCount < obj.getCount())
			{
				highestCount = obj.getCount();
			}
			size++;
		}
		CountObject.steamOutSequenceClose(out);
		out.flush();
		DebugLog.log("Output completed. Unique Rows:"+ size);
		
		Properties props=new Properties();
		props.put("noOfValidRead", totalCounts);
		props.put("lowestCount", lowestCount);
		props.put("highestCount", highestCount);
		
		return props;
	}
	
	private void exchange(int i,int j)
	{
		MinHeapNode temp=arrays[i];
		arrays[i]=arrays[j];
		arrays[j] =temp;
	}
	
//	public static void main(String args[])
//	{
//		//ArraysMergerHeap m=new ArraysMergerHeap("/ifs/scratch/c2b2/hb_lab/dl2868/tmp",158);
//		ArraysMergerHeap m=new ArraysMergerHeap("tmp",21);
//		try
//		{
//			m.output();
//		} catch (Exception e)
//		{
//			// TODO Auto-generated catch block
//			DebugLog.log(e);
//		}
//	}
}
