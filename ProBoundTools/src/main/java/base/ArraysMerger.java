package base;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class ArraysMerger
{
	private int SIZE = 1 * 1000 * 1000; //size might be specified somewhere else
	
	private CountObject[] countObjects;
	private int index = 0;
	private int listIndex = -1;
	private File tempFolder;
	private String filePrefix;
	private int MAX_THREAD_NUM = 2;
	private ExecutorService service = null;
	private Semaphore semaphore = null;
	
	private static Comparator<CountObject> comparator = new SequenceComparator();
	
	public ArraysMerger(String tmpFolder, String filePrefix, int bufferSize) 
	{
		int cores = Util.getNumberOfCores();
		MAX_THREAD_NUM = Util.getMaxThreadNumber();
		
		service = Executors.newFixedThreadPool(MAX_THREAD_NUM);
		semaphore =new Semaphore(MAX_THREAD_NUM+1);
		DebugLog.log("Available cores: "+cores);
		DebugLog.log("MAX_THREAD_NUM:  "+MAX_THREAD_NUM);
		
		try
		{
			//tempFolder = File.createTempFile("selex", "");// SIZE=10;
			//tempFolder = new File("/ifs/scratch/c2b2/hb_lab/dl2868/tmp");
			tempFolder= new File(tmpFolder); //SIZE=10;
			if(bufferSize>0)
			{
				SIZE=bufferSize;
			}
			countObjects = new CountObject[SIZE];
			if (tempFolder.exists())
				tempFolder.delete();
			tempFolder.mkdir();
			
			this.filePrefix = filePrefix;
			
			DebugLog.log("Temp Folder:"+tempFolder);
			DebugLog.log("Temp File Size:"+SIZE);
			
		} 
		catch (Exception e)
		{
			DebugLog.log(e);
		}

		addSpace();
	}

	/**
	 * Exports the CountObject array asynchronously 
	 */
	private void asyncExport()
	{
		final CountObject[] localList = countObjects;
		final int localIndex = listIndex;
		final int size = this.index;
		try
		{
			DebugLog.log("Waiting to submit an export job ... ");
			semaphore.acquire();
		} 
		catch (InterruptedException e)
		{
			DebugLog.log(e);
		}

		// DebugLog.log("asyncExport() done");

		service.submit(new Runnable()
		{

			//@Override
			public void run()
			{
				DebugLog.log("Started for exporting ... " +" Size:"+ (size)+" Index:"+localIndex);
				try
				{
					export(localList, localIndex, size);
				}
				finally
				{
					//DebugLog.log("Releasing semaphore ... ");
					semaphore.release();
				}
				DebugLog.log("Finished exporting ... ");
			}

		});

	}

	/**
	 * Sort the array of CountObjects and output to a binary file
	 * @param localList
	 * @param localIndex
	 * @param listSize
	 */
	private void export(CountObject[] localList, int localIndex, int listSize)
	{
		DebugLog.log("Sorting:");

		sort(localList, listSize);

		DebugLog.log(" Outputing:");
		String filePath = getFileName(tempFolder.getAbsolutePath(), this.filePrefix, localIndex);
		try
		{
			//BufferedOutputStream b=new BufferedOutputStream(new FileOutputStream(filePath));
			
			Long currentTime = System.currentTimeMillis();
			DataOutputStream out = new DataOutputStream(
					new BufferedOutputStream(new FileOutputStream(filePath), 1000*1000*10));
			CountObject[] a = localList;
			int size=0;
			CountObject lastObj = null;
			for (int i = 0; i < listSize; i++)
			{
				CountObject currObj = a[i];
				if(lastObj!=null)
				{
					if(lastObj.compareToByKey(currObj)==0)//collapse the same sequences
					{
						lastObj.setCount(lastObj.getCount()+currObj.getCount());
					}
					else // outputs until there is a different sequence coming up
					{
						CountObject.steamOutSequence(out, lastObj);
						lastObj = currObj;
						size++;
					}
				}
				else
				{
					lastObj = currObj;
				}
			}
			if(lastObj!=null) //outputs the last one
			{
				CountObject.steamOutSequence(out, lastObj);
				size++;
			}
			CountObject.steamOutSequenceClose(out);
			DebugLog.log(" Done output file:" + filePath +" Rows:"+ size);

			out.flush();
			out.close(); // this one is really slow for some reason

			DebugLog.log("File["+filePath+"] closed. Done in ["+((System.currentTimeMillis()- currentTime) / 1000)+"] seconds");
		} catch (Exception e)
		{
			DebugLog.log(e);
			throw new RuntimeException(e);
		}
	}

	private void addSpace()
	{
		if (listIndex != -1)
		{
			asyncExport();
		}
		//DebugLog.log("Refreshing...");
		countObjects=null; //intended for JVM to recycle the space first
		countObjects=new CountObject[SIZE];
		//DebugLog.log("Refreshing...DONE");
		listIndex++;
		index = 0;
	}

	public void add(CountObject obj)
	{
		countObjects[index++] = obj;
		if (index == SIZE)
		{
			//DebugLog.log("Adding more space");
			addSpace();
			//DebugLog.log("Adding more space DONE");
		}
	}

	public void finish()
	{
		asyncExport();
		try
		{
			service.shutdown();
			// now wait for the jobs to finish
			service.awaitTermination(1000000, TimeUnit.SECONDS);

			DebugLog.log("All output done");
		} catch (InterruptedException e)
		{
			// TODO Auto-generated catch block
			DebugLog.log(e);
		}
	}

	/**
	 * Different sorting algorithms(JDK's merge sort, quick sort, radix sort) were tested. 
	 * Turns out Arrays.sort is the fastest most of the time.
	 * @param list
	 * @param size
	 */
	private void sort(CountObject[] list, int size)
	{
		long t1=System.currentTimeMillis();

		/*
		Util.radixSort(list, size);
		*/
		Arrays.sort(list, 0, size ,comparator);
		
		/*
		QuickSort.qsort(list, 0, end, new Comparator<CountObject>()
		{
			@Override
			public int compare(CountObject arg0, CountObject arg1)
			{
				if (arg0 == null)
					return 1;
				if (arg1 == null)
					return -1;
				return arg0.compareToByKey(arg1);
			}

		});*/
		
		long t2=System.currentTimeMillis();
		DebugLog.log("Sorting done in ["+((t2-t1)/1000.0)+"] seconds.");
	}

	/**
	 * Merges sorted files into the final count file.
	 * The sorted files are named with as the combination of prefix + fileNumber
	 * 
	 * @param textOutput
	 * @param path
	 * @return
	 */
	public Properties output(boolean textOutput, String path)
	{
		ArraysMergerHeap heap = new ArraysMergerHeap(
				tempFolder.getAbsolutePath(),this.filePrefix,  this.listIndex + 1);
		try
		{
			Properties props= heap.output(textOutput, path);
			cleanTempFiles();
			return props;
		}
		catch (Exception e)
		{
			// TODO Auto-generated catch block
			DebugLog.log(e);
			return null;
		}
	}
	
	private void cleanTempFiles()
	{
		for(int i=0;i<= this.listIndex;i++)
		{
			 String path= getFileName(tempFolder.getAbsolutePath(), this.filePrefix, i);
			 DebugLog.log("Removing temp file:" + path);
			 new File(path).delete();
		}
	}
	
	public static String getFileName(String folder, String filePrefix, int id)
	{
		return folder + "/" + filePrefix + "-" + id + ".dat"; 
	}
	
}
