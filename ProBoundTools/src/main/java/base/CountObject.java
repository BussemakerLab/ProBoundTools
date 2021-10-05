package base;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class CountObject implements Comparable<CountObject> {

	private Sequence key;
	private int count;
	
    public CountObject(Sequence key, Integer count)
    {
    	this.key=key;
    	this.count=count;
    }

	//@Override
	public int compareTo(CountObject out)
	{
		//DebugLog.log("count:"+count);
		if(this.count==out.count)
		{
			//DebugLog.log("count:"+count+"  "+this.key +" --> "+ out.key +" = " + this.key.compareTo(out.key));
			return this.key.compareTo(out.key);
		}
		return out.count - count ;
	}
	

	public int compareToByKey(CountObject out)
	{
//		if(this.key==null)
//		{
//			return 1;
//		}
//		else if(out==null || out.key==null)
//		{
//			return -1;
//		}
		int c=this.key.compareTo(out.key);
		return c;
	}

	public Sequence getKey()
	{
		return key;
	}

	public void setKey(Sequence key)
	{
		this.key = key;
	}

	public Integer getCount()
	{
		return count;
	}

	public void setCount(Integer count)
	{
		this.count = count;
	}

	@Override
	public String toString()
	{
		return "CountObject [key=" + key + ", count=" + count + "]";
	}
	
	public String toString2()
	{
		return key+ " "+count;
	}
	
	public static void steamOutSequence(DataOutput os, CountObject obj) throws IOException
	{
		//DebugLog.log(obj);
		Sequence seq= (obj.getKey());
		os.writeLong(seq.getValue());
		os.writeInt(seq.getLength());
		os.writeInt((int)obj.getCount());
		//DebugLog.log(seq.getValue()+" "+seq.getLength()+" "+obj.getCount());
	}
	

	public static void steamOutSequenceTxt(DataOutput os, CountObject obj) throws IOException
	{
		//DebugLog.log(obj);
		Sequence seq= (obj.getKey());
		os.write((seq.getString() + "  "+obj.getCount() +"\n").getBytes());
	}
	

	public static void steamOutSequenceClose(DataOutput os) throws IOException
	{
		//DebugLog.log(obj);
		os.writeLong(Long.MAX_VALUE);
		//DebugLog.log(seq.getValue()+" "+seq.getLength()+" "+obj.getCount());
	}
	

	//read in 16 bytes
	public static CountObject steamInSequence(DataInput os) throws IOException
	{
		//DebugLog.log(obj);
		Long v=os.readLong();
		if(v==Long.MAX_VALUE)
		{
			return null;
		}
		Integer l=os.readInt();
		Integer c=os.readInt();
		//DebugLog.log(v+" "+l+" "+c);
		CountObject obj=new CountObject(new Sequence(v,l),c);
		//DebugLog.log("Reading:"+obj);
		return obj;
	}
	
}