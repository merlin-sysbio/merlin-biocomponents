package pt.uminho.ceb.biosystems.merlin.biocomponents.io;

public class ModelSourcesEnumerator {
	
	public enum ModelSources{
		
		MODEL_SEED("ModelSEED smbl model");
		
		private String source;
		
		private ModelSources(String modelSource){
			this.source = modelSource;
		}
		
		public String source(){
			return this.source;
		}
		
	}

}
