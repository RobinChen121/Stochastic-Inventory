package workforce;


public class StaffState {
	int period;
	int iniStaffNum;
		
	public StaffState(int period, int iniStaffNum)
	{
		this.period = period;
		this.iniStaffNum = iniStaffNum;
	}
		
		
	@Override
	public int hashCode(){
		String hash = "";
		hash = hash + period + iniStaffNum;
		return hash.hashCode();
	}
		
	@Override
	public boolean equals(Object o) {
		if (o instanceof StaffState)
			return ((StaffState) o).period == this.period &&
					((StaffState) o).iniStaffNum == this.iniStaffNum;
		else
			return false;
	}
	
	@Override
	public String toString() {
		return "period = " + period +", "+"initial staff number = " + iniStaffNum; 
	}
}
