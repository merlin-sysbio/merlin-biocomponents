package pt.uminho.ceb.biosystems.merlin.biocomponents.io;

public class ModelSourcesEnumerator {
	
	public enum ModelSources{
		
		MODEL_SEED("ModelSEED"),
		
		BIGG("BIGG"){
			@Override
			public String toString(){
				return "BIGG smbl model";
			}
		};
		
		private String source;
		
		private ModelSources(String modelSource){
			this.source = modelSource;
		}
		
		public String sourceName(){
			return this.source;
		}
		
		@Override
		public String toString(){
			return "ModelSEED smbl model";
		}
	}
	
	
	/**
	 * @author amaromorais
	 * available programs for models genomes alignment
	 */
	public static enum AlignmentMethod{	
		
		BLAST,
		SmithWaterman
	}

}
