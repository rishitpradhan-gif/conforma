import {parseProtein} from './proteins.mjs';self.onmessage=e=>{try{self.postMessage({ok:true,structure:parseProtein(e.data.text)});}catch(error){self.postMessage({ok:false,error:error.message});}};
