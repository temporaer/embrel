package src.Loaders {
	
	import flash.events.Event;
	import flash.events.EventDispatcher;
	import flash.events.ProgressEvent;
	import flash.net.URLLoader;
	import flash.net.URLRequest;
	import flash.events.IOErrorEvent;

	
	/**
	* @description XML-Loading-Class for Loading and Providing XML-Data
	* @author Clemens Petzold
	*/
	public class xmlLoader extends EventDispatcher{
		
		public static var XML_LOADED:String = "xmlloaded";
		public static var UPDATE:String = "update";
		public static var LOADINGERROR:String = "onError";
		
		private var _myXML:XML;		
		private var _xmlurl:String = "";
		private var xmlreq:URLRequest;
		private var myLoader:URLLoader;	
		private var _progress:Number = 0;
		
		public function xmlLoader() {
			_myXML = new XML();
			_myXML.ignoreWhitespace = true;
			xmlreq = new URLRequest();
		}		
		
		public function set xmlurl(value:String):void {
			_xmlurl = value;
			xmlreq.url = _xmlurl;
			myLoader = new URLLoader(xmlreq);			
			myLoader.addEventListener("complete", xmlLoaded);
			myLoader.addEventListener(ProgressEvent.PROGRESS, handleProgress);
			myLoader.addEventListener(IOErrorEvent.IO_ERROR, onXMLFileError);
		}
		
		private function handleProgress(e:ProgressEvent):void
		{
			this._progress = e.bytesLoaded / e.bytesTotal;
			
			this.dispatchEvent( new Event( UPDATE ) );
		}
		
		
		private function xmlLoaded(event:Event):void
		{
			_myXML = XML(myLoader.data);
			this.dispatchEvent(new Event(XML_LOADED));			
			myLoader.removeEventListener("complete", xmlLoaded);
			myLoader.removeEventListener(IOErrorEvent.IO_ERROR, onXMLFileError);			
		}		
		
		private function onXMLFileError(e:Event):void
		{
			myLoader.removeEventListener("complete", xmlLoaded);
			myLoader.removeEventListener(IOErrorEvent.IO_ERROR, onXMLFileError);
			trace(this, "XML-Datei " + xmlreq.url + " konnte nicht geladen werden! Stellen Sie sicher das die Datei existiert.");
			this.dispatchEvent( new Event( LOADINGERROR ) );
		}
		
		public function get myXML():XML { return _myXML; }
		
		public function get progress():Number { return _progress; }
		
		private function set progress(value:Number):void {
			_progress = value;
		}
		
	}
	
}