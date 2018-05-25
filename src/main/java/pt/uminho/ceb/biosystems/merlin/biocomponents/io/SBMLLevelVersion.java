/*
 * Copyright 2010
 * IBB-CEB - Institute for Biotechnology and Bioengineering - Centre of Biological Engineering
 * CCTC - Computer Science and Technology Center
 *
 * University of Minho 
 * 
 * This is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * 
 * This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
 * GNU Public License for more details. 
 * 
 * You should have received a copy of the GNU Public License 
 * along with this code. If not, see http://www.gnu.org/licenses/ 
 * 
 * Created inside the SysBioPseg Research Group (http://sysbio.di.uminho.pt)
 */
package pt.uminho.ceb.biosystems.merlin.biocomponents.io;

/**
 * Enumeration of common SBML levels and versions 
 * 
 * @author amaromorais
 */
public enum SBMLLevelVersion {

	L2V1{ // level 2 version 1
		
		public int getVersion(){
			return 1;
		}
		
		public int getLevel(){
			return 2;
		}

		@Override
		public String toString(){
			return "SBML Level 2 version 1";
		}
	}, 
	L2V2{ // level 2 version 2

		public int getVersion(){
			return 2;
		}
		
		public int getLevel(){
			return 2;
		}

		@Override
		public String toString(){
			return "SBML Level 2 version 2";
		}
	},
	L2V3{ //level 2 version 3

		public int getVersion(){
			return 3;
		}
		
		public int getLevel(){
			return 2;
		}

		@Override
		public String toString(){
			return "SBML Level 2 version 3";
		}
	},
	L2V4{ //level 2 version 4 //default

		public int getVersion(){
			return 4;
		}
		
		public int getLevel(){
			return 2;
		}

		@Override
		public String toString(){
			return "SBML Level 2 version 4";
		}
	},
	L3V1{ //level 3 version 1
		
		public int getVersion(){
			return 1;
		}
		
		@Override
		public String toString(){
			return "SBML Level 3 version 1";
		}
	}, 
	L3V2{ //level 3 version 1
		
		public int getVersion(){
			return 2;
		}
		
		@Override
		public String toString(){
			return "SBML Level 3 version 2";
		}
	};

	public int getLevel(){
		return 3;
	}

	public int getVersion(){
		return 2;
	}

	@Override
	public String toString(){
		return "SBML Level 3 version 2";
	}

}
