AFRAME.registerComponent('loading', {
    schema: {
        time: {type: 'number', default: 10},
        show: {type: 'boolean', default: true},
        zDepth: {type: 'number', default: -1}
    },
    init: function() {
        this.initialized = performance.now();
        this.adjusted = false;
        this.width = this.visibleWidthAtZDepth(this.data.zDepth);
        this.height = this.visibleHeightAtZDepth(this.data.zDepth);
        const userAgent = navigator.userAgent || navigator.vendor || window.opera;
        this.isIOS = /iPad|iPhone|iPod/.test(userAgent) && !window.MSStream;
        this.startedPortrait = (matchMedia && matchMedia("(orientation: portrait)").matches) || 
                (["portrait", "portrait-primary", "portrait-secondary"].includes(window.screen.orientation) ||
                (Utils.mobilecheck() && window.height > window.width));
        const geometry = new THREE.PlaneBufferGeometry( this.width, this.height );
        const material = new THREE.MeshBasicMaterial( {color: 0x1F2630, side: THREE.DoubleSide} );
        const plane = new THREE.Mesh( geometry, material );
        this.buffer = document.getElementById("loadingBuffer");
        this.buffer.object3D.add( plane );
    },
    update: function(oldData) {
        if (this.data.show !== oldData.show && !this.data.show) {
            const dismissTime = performance.now();
            const elapsedSeconds = (dismissTime - this.initialized) / 1000;
            if (elapsedSeconds > this.data.time) {
                this.el.object3D.visible = false;
                this.buffer.object3D.visible = false;
                document.getElementById("hud").setAttribute('visible', true);
                document.getElementById("player-camera").setAttribute("look-controls", "");
                document.getElementById("scene").setAttribute("vr-mode-ui", "enabled: true");
            } else {
                setTimeout(() => {
                    this.el.object3D.visible = false;
                    this.buffer.object3D.visible = false;
                    document.getElementById("hud").setAttribute('visible', true);
                    document.getElementById("player-camera").setAttribute("look-controls", "");
                    document.getElementById("scene").setAttribute("vr-mode-ui", "enabled: true");
                }, ((this.data.time - elapsedSeconds) * 1000))
            }
        } 
        const width = this.visibleWidthAtZDepth(this.data.zDepth);
        const height = this.visibleHeightAtZDepth(this.data.zDepth);
        if (!this.adjusted) {
            this.width = width;
            this.height = height;
            this.adjusted = true;
        }
        const matchMedia = window.webkitMatchMedia || window.mozMatchMedia || window.oMatchMedia || window.msMatchMedia || window.matchMedia;
        const isPortrait = (matchMedia && matchMedia("(orientation: portrait)").matches) || 
                            (["portrait", "portrait-primary", "portrait-secondary"].includes(window.screen.orientation) ||
                            (Utils.mobilecheck() && window.height > window.width))
        if (isPortrait) {
            this.el.object3D.rotation.set(0, 0, -32.987);
            if (this.isIOS && this.startedPortrait) {
                this.el.setAttribute('width', this.height);
                this.el.setAttribute('height', this.width);
            } else if (this.isIOS) {
                this.el.setAttribute('width', this.width * .55);
                this.el.setAttribute('height', this.height * .55);
            } else {
                this.el.setAttribute('width', height);
                this.el.setAttribute('height', width);
            }
            document.getElementById('hud').object3D.position.set(-width/2 + .25, height/2 - .25, this.data.zDepth);
        } else {
            this.el.object3D.rotation.set(0, 0, 0);
            if (this.isIOS && this.startedPortrait) {
                this.buffer.setAttribute('geometry', 'width', this.height * 3);
                this.buffer.setAttribute('geometry', 'height', this.width * 3);
                this.buffer.setAttribute('material', 'color', 0x1F2630);
                this.el.setAttribute('width', this.height * 1.5);
                this.el.setAttribute('height', this.width * 1.5);
                console.log(this.el.getAttribute('position').x);
                this.buffer.object3D.position.set(this.el.getAttribute('position').x, this.el.getAttribute('position').y, this.el.getAttribute('position').z - 1);
            } else {
                this.el.setAttribute('width', width);
                this.el.setAttribute('height', height);
            }
            document.getElementById('hud').object3D.position.set(-width/2 + .25, height/2 - .25, this.data.zDepth);
        }
    },
    tick: function() {
        const cameraPosition = document.getElementById('rig').object3D.position;
        this.el.object3D.position.set(cameraPosition.x, cameraPosition.y, cameraPosition.z - 1);
    },
    // These are duplicated, fix later
    visibleHeightAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        // compensate for cameras not positioned at z=0
        const cameraOffset = camera.position.z;
        if ( depth < cameraOffset ) depth -= cameraOffset;
        else depth += cameraOffset;
      
        // vertical fov in radians
        const vFOV = camera.fov * Math.PI / 180;
      
        // Math.abs to ensure the result is always positive
        return 2 * Math.tan( vFOV / 2 ) * Math.abs( depth );
    },
    visibleWidthAtZDepth: function( depth ) {
        const camera = this.el.sceneEl.camera;
        const height = this.visibleHeightAtZDepth( depth, camera );
        let width = height * camera.aspect;
        return width;
    }

});